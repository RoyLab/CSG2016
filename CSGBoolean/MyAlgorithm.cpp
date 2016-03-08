#include <unordered_map>
#include <unordered_set>
#include <list>
#include <queue>

#include "MyAlgorithm.h"
#include "ItstAlg.h"
#include "ItstGraph.h"
#include "LoopletTable.h"
#include "csg.h"
#include "COctree.h"
#include "MyMesh.h"
#include "CGALext.h"

namespace
{
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

    enum Mark { UNVISITED, SEEDED, VISITED };
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

        ItstAlg* itst = new ItstAlg(pMeshList);
        itst->doIntersection(intersectLeaves);
        floodColoring(pCsg, itst);

        SAFE_DELETE(itst);
        SAFE_DELETE(pCsg);
        SAFE_DELETE(pOctree);
        for (MyMesh* pMesh : *pMeshList)
            delete pMesh;
        SAFE_DELETE(pMeshList);
    }

    MyAlgorithm::AutoIndicator MyAlgorithm::computeFullIndicator(VH vh, size_t meshId)
    {
        AutoIndicator ind;
        ind.reset(new Indicator[pMeshList->size()]);

        ind[meshId] = myext::BT_ON;

        if (vh->data && vh->data->proxy)
        {
            const auto& ctx = (*vh->data->proxy).pointer()->ctx;
            for (auto &c : ctx)
                ind[c.meshId] = myext::BT_ON;
        }

        for (size_t meshId = 0; meshId < pMeshList->size(); meshId++)
        {
            if (ind[meshId] == myext::BT_ON) continue;
            ind[meshId] = pointInPolyhedron(vh->point(), pMeshList[meshId]);
        }

    }

    void MyAlgorithm::createFirstSeed(SeedInfoWithMeshId& info)
    {
        VH seedV = (*pMeshList)[info.meshId]->vertices_begin();
        info.indicator = computeFullIndicator(seedV, info.meshId);
        info.seedVertex = seedV;

        info.seedFacet = seedV->halfedge()->facet();
    }

    void MyAlgorithm::floodColoring(CSGTree<MyMesh>* pCsg, ItstAlg* itstAlg)
    {
        typedef CSGTree<MyMesh> CSGTree;

        auto &meshList = *pMeshList;
        size_t nMesh = meshList.size();

        GroupParseInfo infos;

        // 表示该mesh有没有被处理过
        //infos.meshSeedFlag = new bool[nMesh];
        //memset(infos.meshSeedFlag, 0, sizeof(bool));

        for (size_t imesh = 0; imesh < nMesh; imesh++)
        {
            //if (infos.meshSeedFlag[imesh]) continue;
            //infos.meshSeedFlag[imesh] = true;

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

    void MyAlgorithm::floodSimpleGroup(GroupParseInfo& infos)
    {

    }

    void MyAlgorithm::floodComplexGroup(GroupParseInfo& infos, SeedInfo& seed, int* trimTree)
    {
        if (seed.seedFacet->mark == VISITED)
            return;

        int nVIP = trimTree->numberOfVIP();

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
        q.push(s);

        bool first = true;
        std::vector<Loop> loops;

        while (!q.empty())
        {
            if (q.front().fh->mark != VISITED)
            {
                FH fh = q.front().fh;
                fh->mark = VISITED;

                ItstGraph* ig = new ItstGraph;
                ig->createGraph(q.front());
                ig->propInd();

                loops.clear();
                fh->data->iTri->looplets->getAllCircles(loops);

                for (auto& loop : loops)
                {
                    if (ig->classify(loop))
                        addFacet(loop);
                }


                for (size_t i = 0; i < 3; i++)
                {
                    ig->getCornerInds(&s, i);
                    s.fh = ? ;
                    s.fh->mark = SEEDED;

                    if (sameGroup(s.fh, fh))
                        q.push(s);
                    else infos.curMeshSeedQueue.push(s);
                }

                delete ig;
            }

            q.pop();
        }
    }

}