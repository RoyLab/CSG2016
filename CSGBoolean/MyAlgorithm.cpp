#include <unordered_map>
#include <unordered_set>
#include <list>
#include <queue>

#include "MyAlgorithm.h"
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
        pMeshList = new std::vector<MyMesh*>;
        for (MyMesh* pMesh : meshes)
            pMeshList->push_back(new MyMesh(*pMesh));

        CSGTree<MyMesh>* pCsg = new CSGTree<MyMesh>;
        pCsg->createCSGTreeFromExpr(expr, pMeshList->data(), pMeshList->size());
        pCsg->makePositiveAndLeftHeavy();

        Octree *pOctree = new Octree;
        std::vector<Octree::Node*> intersectLeaves;
        pOctree->build(*pMeshList, &intersectLeaves);

        doIntersection(intersectLeaves);
        floodColoring(pCsg);

        SAFE_DELETE(pCsg);
        SAFE_DELETE(pOctree);
        for (MyMesh* pMesh : *pMeshList)
            delete pMesh;
        SAFE_DELETE(pMeshList);
    }

    void MyAlgorithm::doIntersection(std::vector<Octree::Node*>& intersectLeaves)
    {
        typedef std::unordered_set<IndexPair> TriIdSet;
        typedef std::unordered_map<IndexPair, TriIdSet*> MeshIdTriIdMap;

        MeshIdTriIdMap antiOverlapMap;
        antiOverlapMap.max_load_factor(0.6);
        csgResult = new MyMesh;
        meshRelTable = new myext::TriangleTable<bool>(pMeshList->size());
        auto &meshList = *pMeshList;

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

                    // 这里假设map的遍历是保序的，因此meshIdPair自动的分为大小
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

                            meshRelTable->getValue(meshId[0], meshId[1]) = true;

                            setupIsectFacet(fh0);
                            setupIsectFacet(fh1); 

                            if (sign == myext::COPLANAR)
                            {
                                fh0->isectInfo->coplanars.emplace_back(fh1);
                                fh1->isectInfo->coplanars.emplace_back(fh0);
                                continue;
                            }

                            checkNonmanifoldEdge(fh0, fh1, &result, antiOverlapSet);
                            setupPonits(fh0, fh1, result);
                        }
                    }
                }
            }
        }
    }

    MyAlgorithm::AutoIndicator MyAlgorithm::computeFullIndicator(VH fh, size_t meshId)
    {
        AutoIndicator ind;
        ind.reset(new Indicator[pMeshList->size()]);

        ind[meshId] = myext::BT_ON;

        auto &onList = (*fh->shared->agency)->onList;
        for (size_t meshId : onList)
            ind[meshId] = myext::BT_ON;

        for (size_t meshId = 0; meshId < pMeshList->size(); meshId++)
        {
            if (ind[meshId] == myext::BT_ON) continue;
            ind[meshId] = pointInPolyhedron(fh->point(), pMeshList[meshId]);
        }

    }

    void MyAlgorithm::createFirstSeed(SeedInfoWithMeshId& info)
    {
        VH seedV = (*pMeshList)[info.meshId]->vertices_begin();
        info.indicator = computeFullIndicator(seedV, info.meshId);
        info.seedVertex = seedV;

        info.seedFacet = seedV->halfedge()->facet();
    }

    void MyAlgorithm::floodColoring(CSGTree<MyMesh>* pCsg)
    {
        typedef CSGTree<MyMesh> CSGTree;

        auto &meshList = *pMeshList;
        size_t nMesh = meshList.size();

        GroupParseInfo infos;

        // 表示该mesh有没有被处理过
        infos.meshSeedFlag = new bool[nMesh];
        memset(infos.meshSeedFlag, 0, sizeof(bool));

        for (size_t imesh = 0; imesh < nMesh; imesh++)
        {
            if (infos.meshSeedFlag[imesh]) continue;
            infos.meshSeedFlag[imesh] = true;

            SeedInfoWithMeshId seedInfo;
            seedInfo.meshId = imesh;
            createFirstSeed(seedInfo);
            infos.otherMeshSeedQueue.push(seedInfo);

            while (!infos.otherMeshSeedQueue.empty())
            {
                SeedInfoWithMeshId& curSeed = infos.otherMeshSeedQueue.front();
                SeedInfo sndseed;
                sndseed.seedFacet = curSeed.seedFacet;
                sndseed.seedVertex = curSeed.seedVertex;
                CSGTree* fstTrimTree = genFirstTrimTree(curSeed, pCsg, &sndseed.indicator);
                infos.curMeshId = curSeed.meshId;
                infos.curMeshSeedQueue.push(sndseed);

                while (!infos.curMeshSeedQueue.empty())
                {
                    SeedInfo& sndInfo = infos.curMeshSeedQueue.front();
                    AutoIndicator ind;
                    CSGTree* sndTrimTree = genSecondTrimTree(sndInfo, fstTrimTree, &ind);

                    if (!ind) floodSimpleGroup(sndTrimTree, sndInfo, infos);
                    else floodComplexGroup(sndTrimTree, sndInfo, ind, infos);

                    infos.curMeshSeedQueue.pop();
                }

                infos.otherMeshSeedQueue.pop();
            }
        }
    }

    void MyAlgorithm::floodSimpleGroup()
    {

    }

    void MyAlgorithm::floodComplexGroup()
    {
        typedef CSGTree<MyMesh> CSGTree;
        GroupParseInfo infos;
        SeedInfo seed;
        CSGTree* trimTree;

        if (seed.seedFacet->mark == VISITED) return;

        size_t nVIP = trimTree->numberOfVIP();

        struct Seed
        {
            FH fh;
            VH vh;
            Indicator *ind = nullptr;
        } s;

        Queue<Seed> q;         
        s.fh = seed.seedFacet;
        s.vh = seed.seedVertex;
        s.ind = new Indicator[nVIP];



        q.push(seed.seedFacet);

        while (!q.empty())
        {
            if (q.front()->mark != VISITED)
            {
                q.front()->mark = VISITED;
                createGraph(q.front());

            }

            q.pop();
        }
    }

}