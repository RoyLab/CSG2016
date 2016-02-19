#include <unordered_map>
#include <unordered_set>
#include <list>

#include "MyAlgorithm.h"
#include "SharedZone.h"
#include "csg.h"
#include "COctree.h"
#include "MyMesh.h"
#include "CGALext.h"

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

    template <class Itr>
    inline bool operator>(Itr& a, Itr& b)
    {
        return &(*a) > &(*b);
    }

    template <class Itr>
    inline bool operator<(Itr& a, Itr& b)
    {
        return &(*a) < &(*b);
    }

    template <class Itr>
    inline bool operator>=(Itr& a, Itr& b)
    {
        return &(*a) >= &(*b);
    }

    template <class Itr>
    inline bool operator<=(Itr& a, Itr& b)
    {
        return &(*a) <= &(*b);
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
        //szone = new SharedZone(meshList);
        csgResult = new MyMesh;

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

                    TriIdSet* antiOverlapSet = nullptr;
                    auto searchRes = antiOverlapMap.find(meshIdPair);
                    if (searchRes == antiOverlapMap.end())
                    {
                        antiOverlapSet = new TriIdSet;
                        antiOverlapSet->max_load_factor(0.6);
                        antiOverlapMap.emplace(meshIdPair, antiOverlapSet);
                    }
                    else antiOverlapSet = searchRes->second;

                    for (MyMesh::Face_handle fh0 : *meshItr[0]->second)
                    {
                        for (MyMesh::Face_handle fh1 : *meshItr[1]->second)
                        {
                            uint32_t triId[2] = { fh0->id(), fh1->id() };
                            IndexPair triIdPair;
                            MakeIndex(triId, triIdPair);

                            if (antiOverlapSet->find(triIdPair) != antiOverlapSet->end()) 
                                continue;
                            antiOverlapSet->insert(triIdPair);

                            myext::TriTriIsectResult<K> result;
                            myext::Sign sign = myext::tri_tri_intersect(fh0->triangle, fh0->normal,
                                fh1->triangle, fh1->normal, &result);

                            if (sign == myext::NOT_INTERSECT || sign == myext::INTERSECT_ON_POINT) 
                                continue;

                            setupIsectFacet(fh0);
                            setupIsectFacet(fh1); 

                            if (sign == myext::COPLANAR)
                            {
                                fh0->isectInfo->coplanars.emplace_back(fh1);
                                fh1->isectInfo->coplanars.emplace_back(fh0);
                                continue;
                            }

                            // 检查是否是non-manifold的情况，去重
                            checkNonmanifoldEdge(fh0, fh1, &result, antiOverlapSet);
                             
                            setupPonits(fh0, fh1, result);
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