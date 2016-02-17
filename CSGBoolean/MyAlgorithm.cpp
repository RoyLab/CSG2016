#include "MyAlgorithm.h"
#include "csg.h"
#include "COctree.h"
#include "MyMesh.h"

#include <unordered_map>
#include <unordered_set>


namespace
{
    typedef uint64_t IndexPair;

    static inline void MakeIndex(const uint32_t id[], IndexPair& indexPair)
    {
        indexPair = id[1];
        indexPair = indexPair << 32;
        indexPair |= id[0];
    }

    static inline void GetIDFromIndex(uint32_t ID[], const IndexPair& indexPair)
    {
        ID[0] = indexPair & 0xffffffff;
        ID[1] = indexPair >> 32;
    }
}

namespace CSG
{
    void MyAlgorithm::solve(const std::string& expr, std::vector<MyMesh*>& meshes)
    {
        std::vector<MyMesh*> meshList;
        for (MyMesh* pMesh : meshes)
            meshList.push_back(new MyMesh(*pMesh));

        CSGTree<MyMesh>* pCsg = new CSGTree<MyMesh>;
        pCsg->createCSGTreeFromExpr(expr, meshList.data(), meshList.size());
        pCsg->makePositiveAndLeftHeavy();

        Octree *pOctree = new Octree;
        std::vector<Octree::Node*> intersectLeaves;
        pOctree->build(meshList, &intersectLeaves);

        doIntersection(meshList, intersectLeaves);
        floodColoring(pCsg, pOctree);

        SAFE_DELETE(pCsg);
        SAFE_DELETE(pOctree);
        for (MyMesh* pMesh : meshList)
            delete pMesh;
    }


    void MyAlgorithm::doIntersection(std::vector<MyMesh*>& meshList, std::vector<Octree::Node*>& intersectLeaves)
    {
        typedef std::unordered_set<IndexPair> TriIdSet;
        typedef std::unordered_map<IndexPair, TriIdSet*> MeshIdTriIdMap;

        MeshIdTriIdMap antiOverlapMap;
        antiOverlapMap.max_load_factor(0.6);

        for (Octree::Node* leaf : intersectLeaves)
        {
            auto iEnd = leaf->triTable.cend();
            decltype(leaf->triTable.begin()) meshItr[2];
            for (meshItr[0] = leaf->triTable.begin(); meshItr[0] != iEnd; ++meshItr[0])
            {
                meshItr[1] = meshItr[0]; ++meshItr[1];
                for (; meshItr[1] != iEnd; ++meshItr[1])
                {
                    uint32_t meshId[2] = { meshItr[0]->first, meshItr[1]->first };
                    MyMesh* meshes[2] = { meshList[meshId[0]], meshList[meshId[1]] };

                    // 这里假设map的遍历是保序的
                    IndexPair meshIdPair;
                    MakeIndex(meshId, meshIdPair);
                    auto &antiOverlapSet = antiOverlapMap[meshIdPair];

                    for (MyMesh::Face_handle fh0 : *meshItr[0]->second)
                    {
                        for (MyMesh::Face_handle fh1 : *meshItr[1]->second)
                        {
                            uint32_t triId[2] = { fh0->id(), fh1->id() };
                            IndexPair triIdPair;
                            MakeIndex(triId, triIdPair);

                            if (antiOverlapSet->find(triIdPair) != antiOverlapSet->end()) continue;
                            else antiOverlapSet->insert(triIdPair);


                            //auto itr = leaf->TriangleTable.begin();
                            //decltype(iEnd) itr2;
                            //unsigned i, j, ni, nj;
                            //MPMesh *meshi, *meshj;
                            //MPMesh::FaceHandle tri1, tri2;
                            //Vec3d *v0, *v1, *v2, nv, *u0, *u1, *u2, nu, start, end;
                            //MPMesh::FVIter fvItr;
                            //int isISect;
                            //int startT(0), endT(0);
                            //ISectTriangle **si = nullptr, **sj = nullptr;
                            //VertexPos startiT, startjT, endiT, endjT;
                            //ISVertexItr vP1, vP2;

                            // intersection test main body
                            GetCorners(meshi, tri1, v0, v1, v2);
                            GetCorners(meshj, tri2, u0, u1, u2);

                            nv = meshi->normal(tri1);
                            nu = meshj->normal(tri2);

                            startT = INNER; endT = INNER; // return to Zero.

                            //auto &isecTris = (*si)->isecTris[meshj->ID];
                            //auto &icoplTris = (*si)->coplanarTris;

                            isISect = TriTriIntersectTest(*v0, *v1, *v2, nv,
                                *u0, *u1, *u2, nu, startT, endT, start, end);

                            if (isISect < 0) continue;

                            si = &meshi->property(meshi->SurfacePropHandle, tri1);
                            sj = &meshj->property(meshj->SurfacePropHandle, tri2);

                            if (!*si) *si = new ISectTriangle(meshi, tri1);
                            if (!*sj) *sj = new ISectTriangle(meshj, tri2);

                            if (isISect == 0)
                            {
                                (*si)->coplanarTris[meshj->ID].emplace_back(tri2);
                                (*sj)->coplanarTris[meshi->ID].emplace_back(tri1);
                                continue;
                            }

                            startiT = VertexPos(startT & 0xffff);
                            startjT = VertexPos(startT >> 16);

                            endiT = VertexPos(endT & 0xffff);
                            endjT = VertexPos(endT >> 16);

                            if (IsEqual(start, end))
                            {
                                // 点相交
                                Vec3d point = (start + end) / 2;

                                // 最后一个参数表示，可能存在两个以上的插入点
                                vP1 = InsertPoint(*si, startiT, point);
                                InsertPoint(*si, endiT, vP1);
                                InsertPoint(*sj, startjT, vP1);
                                InsertPoint(*sj, endjT, vP1);
                            }
                            else
                            {
                                // 线相交
                                double d = OpenMesh::dot(OpenMesh::cross(nv, nu), end - start);
                                if (meshi->bInverse ^ meshj->bInverse) d = -d;
                                if (d > 0)
                                {
                                    vP1 = InsertPoint(*si, startiT, start);
                                    vP2 = InsertPoint(*si, endiT, end);
                                    InsertSegment(*si, vP1, vP2, *sj);
                                    vP1 = InsertPoint(*sj, startjT, vP1);
                                    vP2 = InsertPoint(*sj, endjT, vP2);
                                    InsertSegment(*sj, vP2, vP1, *si);
                                }
                                else
                                {
                                    vP1 = InsertPoint(*si, startiT, start);
                                    vP2 = InsertPoint(*si, endiT, end);
                                    InsertSegment(*si, vP2, vP1, *sj);
                                    vP1 = InsertPoint(*sj, startjT, vP1);
                                    vP2 = InsertPoint(*sj, endjT, vP2);
                                    InsertSegment(*sj, vP1, vP2, *si);
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    void MyAlgorithm::floodColoring(CSGTree<MyMesh>* pCsg, Octree* pOctree)
    {

    }
}