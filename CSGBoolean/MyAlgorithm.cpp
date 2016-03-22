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
#include "Octree.h"
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

    enum Mark { SEEDED, VISITED, UNVISITED = MARK_BEGIN };
}

namespace CSG
{
    void MyAlgorithm::solve(const std::string& expr, std::vector<MyMesh*>& meshes)
    {
        pMeshList = &meshes;

        CSGTree<MyMesh>* pCsg = new CSGTree<MyMesh>;
        pCsg->createCSGTreeFromExpr(expr, pMeshList->data(), pMeshList->size());
        pCsg->makePositiveAndLeftHeavy();

        pOctree = new Octree;
        std::vector<Octree::Node*> intersectLeaves;
        pOctree->build(*pMeshList, &intersectLeaves);

        ItstAlg* itst = new ItstAlg(pMeshList);
        itst->doIntersection(intersectLeaves);

#ifdef _DEBUG
        itst->computeDebugInfo();
#endif
        floodColoring(pCsg, itst);

        SAFE_DELETE(itst);
        SAFE_DELETE(pCsg);
        SAFE_DELETE(pOctree);
    }

    IIndicatorVector* MyAlgorithm::computeFullIndicator(VH vh, size_t meshId)
    {
        IIndicatorVector *pind = new FullIndicatorVector(pMeshList->size());
        IIndicatorVector &ind = *pind;
        ind[meshId] = REL_SAME;

        if (vh->data && vh->data->proxy)
        {
            const auto& ctx = (*vh->data->proxy).pointer()->ctx;
            for (auto &c : ctx)
                ind[c.meshId] = REL_SAME;
        }

        for (size_t meshId = 0; meshId < pMeshList->size(); meshId++)
        {
            if (ind[meshId] == REL_SAME)
                continue;
            else
                ind[meshId] = PolyhedralInclusionTest(vh->point(), pOctree, *pMeshList, meshId, false);
        }
        return pind;

    }

    void MyAlgorithm::createFirstSeed(SeedInfoWithId& info)
    {
        VH seedV = (*pMeshList)[info.meshId]->vertices_begin();

        info.indicators.reset(computeFullIndicator(seedV, info.meshId));
        info.seedVertex = seedV;
        info.seedFacet = seedV->halfedge()->facet();
        assert(seedV->halfedge()->vertex() == seedV);
    }

    void MyAlgorithm::floodColoring(CSGTree<MyMesh>* pCsg, ItstAlg* itstAlg)
    {
        auto &meshList = *pMeshList;
        size_t nMesh = meshList.size();

        GroupParseInfo infos;
        std::vector<int> itstPrims;

        for (size_t imesh = 0; imesh < nMesh; imesh++)
        {
            SeedInfoWithId idSeed;
            idSeed.meshId = imesh;
            createFirstSeed(idSeed);
            infos.otherMeshSeedQueue.push(idSeed);

            while (!infos.otherMeshSeedQueue.empty())
            {
                SeedInfoWithId& curSeed = infos.otherMeshSeedQueue.front();
                itstPrims.clear();
                itstAlg->get_adjGraph()->getIntersectPrimitives(curSeed.meshId, itstPrims);// unimplemented

            //    //infos.ttree1.reset(new TrimCSGTree<MyMesh>(*pCsg, *curSeed.indicators, itstPrims));
                infos.curMeshId = curSeed.meshId;

            //    hintSeed.indicators = infos.ttree1->downcast(*curSeed.indicators);
                infos.curMeshSeedQueue.push(curSeed);

                while (!infos.curMeshSeedQueue.empty())
                {
                    SeedInfo& sndInfo = infos.curMeshSeedQueue.front();
                    if (sndInfo.seedFacet->mark != VISITED)
                    {
                        //infos.ttree2.reset(new TrimCSGTree<MyMesh>(*infos.ttree1, *sndInfo.indicators, sndInfo.seedFacet));
                        if (!sndInfo.seedFacet->isSimple()) 
                            floodSimpleGroup(infos, sndInfo);
                        else 
                            floodComplexGroup(infos, sndInfo);

                        sndInfo.seedFacet->mark = VISITED;
                        std::cout << 1;
                    }
                    infos.curMeshSeedQueue.pop();
                }

                infos.otherMeshSeedQueue.pop();
            }
        }
    }

    Relation MyAlgorithm::relationOfContext(FH coins, VH vh)
    {
        return REL_NOT_AVAILABLE;
    }

    Relation MyAlgorithm::relationOfContext(Context<MyMesh>& ctx, VH vh, FH &coins)
    {
        CGAL::Oriented_side side;
        switch (ctx.type)
        {
        case CT_VERTEX:
            // TODO
        case CT_EDGE:
            // TODO
        case CT_FACET:
            side = (*ctx.fh)->data->sp.oriented_side(vh->point());
            switch (side)
            {
            case CGAL::ON_NEGATIVE_SIDE:
                return REL_INSIDE;
            case CGAL::ON_POSITIVE_SIDE:
                return REL_OUTSIDE;
            case CGAL::ON_ORIENTED_BOUNDARY:
                coins = *ctx.fh;
                return REL_SAME;
            default:
                ReportError();
                break;
            }
            break;
        default:
            ReportError();
            break;
        }
        return REL_NOT_AVAILABLE;
    }

    void MyAlgorithm::figureOutFaceInds(VH p, VH q, VH r, int meshId, IIndicatorVector* inds)
    {
        auto &context = p->data->proxy->pointer()->ctx;
        FH coincident;
        assert(!CGAL::collinear(p->point(), q->point(), r->point()));

        for (auto &ctx : context)
        {
            if (ctx.meshId == meshId) continue;
            Relation rel = relationOfContext(ctx, q, coincident);
            switch (rel)
            {
            case CSG::REL_INSIDE:
            case CSG::REL_OUTSIDE:
                break;
            case CSG::REL_SAME:
                rel = relationOfContext(coincident, r);
                break;
            default:
                ReportError();
                break;
            }

            if (pMeshList->at(ctx.meshId)->bInverse)
                inverseRelation(rel);

            (*inds)[ctx.meshId] = rel;
        }
    }

    void MyAlgorithm::figureOutFaceInds(SeedInfo& s, int meshId)
    {
        VH pt2 = s.seedVertex->halfedge()->next()->vertex();
        VH pt3 = pt2->halfedge()->next()->vertex();

        figureOutFaceInds(s.seedVertex, pt2, pt3, meshId, s.indicators.get());
    }

    void MyAlgorithm::floodSimpleGroup(GroupParseInfo& infos, SeedInfo& s)
    {
        // 找到这一群的共同indicators
        if (s.seedVertex->data && s.seedVertex->data->proxy->pointer()->ctx.size() > 1)
            figureOutFaceInds(s, infos.curMeshId);

        Queue<FH> queue;
        queue.push(s.seedFacet);

    //    bool isOn = infos.ttree1.eval(s.indicators);

    //    while (!queue.empty())
    //    {
    //        if (queue.front()->mark != VISITED)
    //        {
    //            FH curface = queue.front();
    //            curface->mark = VISITED;
    //            if (isOn) AddFacet(curface);

    //            auto itr = curface->facet_begin();
    //            for (int i = 0; i < 3; i++)
    //            {
    //                if (itr->facet()->mark < UNVISITED) // 跟ray-tracing有关
    //                    continue;

    //                if (isSameGroup(itr->facet(), curface))
    //                    queue.push(itr->facet());
    //                else
    //                {
    //                    SeedInfoWithHint seed2;
    //                    seed2.indicators = infos.ttree1->createIndicatorVector();
    //                    seed2.seedVertex = itr->vertex();
    //                    seed2.seedFacet = itr->facet();

    //                    s.indicators->copy(*seed2.indicators);
    //                    includeOnCase(seed2);
    //                    createHint(seed2);

    //                    infos.curMeshSeedQueue.push(seed2);
    //                }

    //                itr->facet()->mark == SEEDED;
    //                itr = itr->next;
    //            }
    //        }
    //        queue.pop();
    //    }
    }

    void MyAlgorithm::floodComplexGroup(GroupParseInfo& infos, SeedInfo& s)
    {
        //    int nVIP = trimTree->numberOfVIP();

        //    SeedInfo s;

        //    Queue<SeedInfo> q;
        //    s.fh = seed.seedFacet;
        //    s.vh = seed.seedVertex;
        //    s.ind = new Indicator[nVIP];
        //    q.push(s);

        //    bool first = true;
        //    std::vector<Loop> loops;

        //    while (!q.empty())
        //    {
        //        if (q.front().fh->mark != VISITED)
        //        {
        //            FH fh = q.front().fh;
        //            fh->mark = VISITED;

        //            ItstGraph* ig = new ItstGraph;
        //            ig->createGraph(q.front());
        //            ig->propInd();

        //            loops.clear();
        //            fh->data->iTri->looplets->getAllCircles(loops);

        //            for (auto& loop : loops)
        //            {
        //                if (ig->classify(loop))
        //                    addFacet(loop);
        //            }

        //            for (size_t i = 0; i < 3; i++)
        //            {
        //                ig->getCornerInds(&s, i);
        //                s.fh = ? ;
        //                s.fh->mark = SEEDED;

        //                if (IsSameGroup(s.fh, fh))
        //                    q.push(s);
        //                else
        //                {
        //                    
        //                }
        //            }

        //            delete ig;
        //        }

        //        q.pop();
        //    }
        //}
    }

}