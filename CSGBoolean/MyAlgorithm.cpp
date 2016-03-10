#include <unordered_map>
#include <unordered_set>
#include <list>
#include <queue>

#include "MyAlgorithm.h"
#include "ItstAlg.h"
#include "ItstGraph.h"
#include "LoopletTable.h"
#include "TrimCSGTree.h"
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

    IndicatorVector* MyAlgorithm::computeFullIndicator(VH vh, size_t meshId)
    {
        IndicatorVector *pind = new IndicatorVector(pMeshList->size());
        IndicatorVector &ind = *pind;
        ind[meshId] = myext::BT_ON;

        if (vh->data && vh->data->proxy)
        {
            const auto& ctx = (*vh->data->proxy).pointer()->ctx;
            for (auto &c : ctx)
                ind[c.meshId] = myext::BT_ON;
        }

        for (size_t meshId = 0; meshId < pMeshList->size(); meshId++)
        {
            if (ind[meshId] == myext::BT_ON)
                continue;
            else
                ind[meshId] = pointInPolyhedron(vh->point(), pMeshList[meshId]);
        }
        return pind;

    }

    void MyAlgorithm::createFirstSeed(SeedInfo& info)
    {
        VH seedV = (*pMeshList)[info.meshId]->vertices_begin();

        info.indicators = computeFullIndicator(seedV, info.meshId);
        info.seedVertex = seedV;
        info.seedFacet = seedV->halfedge()->facet();
    }

    void MyAlgorithm::floodColoring(CSGTree<MyMesh>* pCsg, ItstAlg* itstAlg)
    {
        auto &meshList = *pMeshList;
        size_t nMesh = meshList.size();

        GroupParseInfo infos;
        std::vector<int> itstPrims;

        // 表示该mesh有没有被处理过
        //infos.meshSeedFlag = new bool[nMesh];
        //memset(infos.meshSeedFlag, 0, sizeof(bool));

        for (size_t imesh = 0; imesh < nMesh; imesh++)
        {
            //if (infos.meshSeedFlag[imesh]) continue;
            //infos.meshSeedFlag[imesh] = true;

            SeedInfoWithId idSeed;
            idSeed.meshId = imesh;
            createFirstSeed(idSeed);
            infos.otherMeshSeedQueue.push(idSeed);

            while (!infos.otherMeshSeedQueue.empty())
            {
                SeedInfoWithId& curSeed = infos.otherMeshSeedQueue.front();
                itstPrims.clear();
                itstAlg->get_adjGraph()->getIntersectPrimitives(curSeed.meshId, itstPrims);

                infos.ttree1.reset(new TrimCSGTree<MyMesh>(*pCsg, *curSeed.indicators, itstPrims));
                infos.curMeshId = curSeed.meshId;

                SeedInfoWithHint hintSeed;
                hintSeed.seedFacet = curSeed.seedFacet;
                hintSeed.seedVertex = curSeed.seedVertex;
                hintSeed.indicators = infos.ttree1->downcast(*curSeed.indicators);

                infos.curMeshSeedQueue.push(hintSeed);

                while (!infos.curMeshSeedQueue.empty())
                {
                    SeedInfoWithHint& sndInfo = infos.curMeshSeedQueue.front();
                    if (sndInfo.seedFacet->mark != VISITED)
                    {
                        infos.ttree2.reset(new TrimCSGTree<MyMesh>(*infos.ttree1, *sndInfo.indicators, sndInfo.seedFacet));
                        if (!sndInfo.seedFacet->isSimple()) floodSimpleGroup(infos, sndInfo);
                        else floodComplexGroup(infos, sndInfo);
                    }
                    infos.curMeshSeedQueue.pop();
                }

                infos.otherMeshSeedQueue.pop();
            }
        }
    }

    void MyAlgorithm::floodSimpleGroup(GroupParseInfo& infos, SeedInfoWithHint& s)
    {
        if (s.seedVertex->isShared())
            correctSeedIndVec(s);

        Queue<FH> queue;
        queue.push(s.seedFacet);

        bool isOn = infos.ttree1.eval(s.indicators);

        while (!queue.empty())
        {
            if (queue.front()->mark != VISITED)
            {
                FH curface = queue.front();
                curface->mark = VISITED;
                if (isOn) AddFacet(curface);

                auto itr = curface->facet_begin();
                for (int i = 0; i < 3; i++)
                {
                    if (itr->facet()->mark != UNVISITED)
                        continue;

                    if (isSameGroup(itr->facet(), curface))
                        queue.push(itr->facet());
                    else
                    {
                        SeedInfoWithHint seed2;
                        seed2.indicators = infos.ttree1->createIndicatorVector();
                        seed2.seedVertex = itr->vertex();
                        seed2.seedFacet = itr->facet();

                        s.indicators->copy(*seed2.indicators);
                        correctSeedIndVec(seed2.indicators, seed2.seedVertex);

                        createHint(seed2);

                        infos.curMeshSeedQueue.push(seed2);
                    }

                    itr->facet()->mark == SEEDED;
                    itr = itr->next;
                }
            }
            queue.pop();
        }
    }

    void MyAlgorithm::floodComplexGroup(GroupParseInfo& infos, SeedInfo& seed, int* trimTree)
    {
        int nVIP = trimTree->numberOfVIP();

        SeedInfo s;

        Queue<SeedInfo> q;
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

                    if (IsSameGroup(s.fh, fh))
                        q.push(s);
                    else
                    {
                        
                    }
                }

                delete ig;
            }

            q.pop();
        }
    }

}