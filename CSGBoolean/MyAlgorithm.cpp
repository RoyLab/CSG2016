#include <unordered_map>
#include <unordered_set>
#include <list>
#include <queue>
#include <deque>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include "MyAlgorithm.h"
#include "ItstAlg.h"
#include "ItstGraph.h"
#include "LoopletTable.h"
#include "TrimCSGTree.h"
#include "csg.h"
#include "Octree.h"
#include "MyMesh.h"
#include "CGALext.h"
#include "BinaryTree.h"
#include "BSPTree.h"

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

        for (int i = 0; i < 8; i++)
            seedV++;

        info.indicators.reset(computeFullIndicator(seedV, info.meshId));
        info.seedVertex = seedV;
        info.seedFacet = seedV->halfedge()->facet();
        assert(seedV->halfedge()->vertex() == seedV);
        assert(info.seedFacet == info.seedVertex->halfedge()->facet());

        // 这条并不成立
        //assert(info.seedVertex == info.seedFacet->halfedge()->vertex());
    }

    void MyAlgorithm::floodColoring(CSGTree<MyMesh>* pCsg, ItstAlg* itstAlg)
    {
        auto &meshList = *pMeshList;
        size_t nMesh = meshList.size();
        csgResult.reset(new MyMesh);
        pConstruct = new Delegate<HalfedgeDS>;

        GroupParseInfo infos;
        std::vector<int> itstPrims;

        auto tree = pCsg->auxiliary();
        infos.curTreeLeaves = new CSGTreeNode*[nMesh];

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
                itstAlg->get_adjGraph()->getIntersectPrimitives(curSeed.meshId, itstPrims);

                infos.ttree1.reset(copy2(tree->pRoot, infos.curTreeLeaves));
                infos.curMeshId = curSeed.meshId;
                infos.curMeshSeedQueue.push(curSeed);

                while (!infos.curMeshSeedQueue.empty())
                {
                    SeedInfo& sndInfo = infos.curMeshSeedQueue.front();
                    if (sndInfo.seedFacet->mark != VISITED)
                    {
                        if (sndInfo.seedFacet->isSimple()) 
                            floodSimpleGroup(infos, sndInfo);
                        else 
                            floodComplexGroup(infos, sndInfo);

                        sndInfo.seedFacet->mark = VISITED;
                    }
                    infos.curMeshSeedQueue.pop();
                }
                infos.otherMeshSeedQueue.pop();
            }
        }
        csgResult->delegate(*pConstruct);
        SAFE_DELETE(pConstruct);
        SAFE_DELETE_ARRAY(infos.curTreeLeaves);
    }

    Relation MyAlgorithm::relationOfContext(FH coins, VH vh)
    {
        return REL_NOT_AVAILABLE;
    }

    Relation convert(CGAL::Oriented_side side)
    {
        switch (side)
        {
        case CGAL::ON_NEGATIVE_SIDE:
            return REL_INSIDE;
        case CGAL::ON_POSITIVE_SIDE:
            return REL_OUTSIDE;
        case CGAL::ON_ORIENTED_BOUNDARY:
            return REL_SAME;
        default:
            ReportError();
            return REL_NOT_AVAILABLE;
        }
    }

    Relation determineEdgeBSP(MyMesh::Halfedge_handle eh, MyMesh::Vertex_handle vh)
    {
        auto f0 = eh->facet();
        auto f1 = eh->opposite()->facet();
        
        auto r0 = f0->data->sp.oriented_side(vh->point());
        auto r1 = f1->data->sp.oriented_side(vh->point());
        
        if (r0 == r1) return convert(r0);
        
        if (r0 == REL_SAME)
        {
            auto opvh = eh->next()->vertex();
            auto r01 = f1->data->sp.oriented_side(opvh->point());
            if (r01 == r1)
                return convert(r0);
            else
                return convert(r1);
        }
        else
        {
            auto opvh = eh->opposite()->next()->vertex();
            auto r10 = f0->data->sp.oriented_side(opvh->point());
            if (r10 == r0)
                return convert(r1);
            else
                return convert(r0);
        }
    }

    Relation determineVertexBSP(MyMesh::Vertex_handle ctx, MyMesh::Vertex_handle vh)
    {
        typedef MyMesh::Face_handle FH;
        std::vector<FH> facets;
        auto end = ctx->halfedge();
        auto cur = end;
        do
        {
            facets.push_back(cur->facet());
            cur = cur->next()->opposite();
        } while (cur != end);

        BSPTree* pTree = new BSPTree(facets);
        SAFE_DELETE(pTree);

        return REL_NOT_AVAILABLE;
    }

    Relation MyAlgorithm::relationOfContext(Context<MyMesh>& ctx, VH vh, FH &coins)
    {
        CGAL::Oriented_side side, side1;
        Relation result = REL_NOT_AVAILABLE;
        switch (ctx.type)
        {
        case CT_VERTEX:
            result = determineVertexBSP(*ctx.vh, vh);
        case CT_EDGE:
            result = determineEdgeBSP(*ctx.eh, vh);
            break;
        case CT_FACET:
            side = (*ctx.fh)->data->sp.oriented_side(vh->point());
            result =  convert(side);
            if (result == REL_SAME) coins = (*ctx.fh);
            break;
        default:
            ReportError();
            break;
        }
        return result;
    }

    void MyAlgorithm::figureOutFaceInds(VH p, VH q, VH r, int meshId, IIndicatorVector* pinds)
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

            (*pinds)[ctx.meshId] = rel;
        }
    }

    void MyAlgorithm::figureOutFaceInds(SeedInfo& s, int meshId)
    {
        ReportError("Not fully implement");
        VH pt2 = s.seedVertex->halfedge()->next()->vertex();
        VH pt3 = pt2->halfedge()->next()->vertex();

        figureOutFaceInds(s.seedVertex, pt2, pt3, meshId, s.indicators.get());
    }

    void MyAlgorithm::AddFacet(FH fh)
    {
        for (int i = 0; i < 3; i++)
        {
            if (fh->vertices[i]->idx == -1)
            {
                if (fh->vertices[i]->data && fh->vertices[i]->data->proxy->pointer()->idx != -1)
                    fh->vertices[i]->idx = fh->vertices[i]->data->proxy->pointer()->idx;
                else
                {
                    int idx = pConstruct->add_vertex(fh->vertices[i]->point());
                    fh->vertices[i]->idx = idx;
                    if (fh->vertices[i]->data)
                        fh->vertices[i]->data->proxy->pointer()->idx = idx;
                }
            }
        }

        pConstruct->addFacets(fh->vertices[0]->idx,
            fh->vertices[1]->idx, fh->vertices[2]->idx);
    }

    bool MyAlgorithm::isSameGroup(FH fh0, FH fh1) const
    {
        if (fh0->isSimple() ^ fh1->isSimple())
            return false;

        if (fh0->isSimple()) return true;
        return fh0->data->itstTri->meshIds == fh1->data->itstTri->meshIds;
    }

    void MyAlgorithm::genVertexInds(IIndicatorVector* target, VH vh) const
    {
        auto &t = *target;
        if (!vh->data || vh->data->proxy->pointer()->ctx.size() == 0)
            return;

        for (auto& ctx : vh->data->proxy->pointer()->ctx)
        {
            t[ctx.meshId] = REL_SAME;
        }
    }

    void MyAlgorithm::floodSimpleGroup(GroupParseInfo& infos, SeedInfo& s)
    {
        // 找到这一群的共同indicators
        if (s.seedVertex->data && s.seedVertex->data->proxy->pointer()->ctx.size() > 1)
            figureOutFaceInds(s, infos.curMeshId);

        Queue<FH> queue;
        queue.push(s.seedFacet);

        Relation *curRelationTable = new Relation[pMeshList->size()];
        for (size_t i = 0; i < pMeshList->size(); i++)
            curRelationTable[i] = static_cast<Relation>((*s.indicators)[i]);

        TestTree testList;
        CSGTreeNode* curTree = copy2(infos.ttree1.get(), infos.curTreeLeaves);
        auto curRelation = ParsingCSGTree((*pMeshList)[infos.curMeshId], curRelationTable,
            pMeshList->size(), curTree, infos.curTreeLeaves, testList);
        assert(!testList.size() || !testList.begin()->testTree->Parent);

        bool needAdd = false;
        if (curRelation == REL_SAME)
            needAdd = true;

        needAdd = true;

        while (!queue.empty())
        {
            if (queue.front()->mark != VISITED)
            {
                FH curface = queue.front();
                curface->mark = VISITED;
                if (needAdd) AddFacet(curface);

                for (int i = 0; i < 3; i++)
                {
                    auto neighbour = curface->edges[i]->opposite()->facet();
                    if (neighbour->mark == SEEDED || neighbour->mark == VISITED) // 跟ray-tracing有关
                        continue;

                    if (isSameGroup(neighbour, curface))
                        queue.push(neighbour);
                    else
                    {
                        SeedInfo seed2;
                        seed2.seedFacet = neighbour;
                        seed2.seedVertex = curface->edges[i]->opposite()->vertex();
                        seed2.indicators.reset(new FullIndicatorVector(
                            *reinterpret_cast<FullIndicatorVector*>(s.indicators.get())));
                        genVertexInds(seed2.indicators.get(), seed2.seedVertex);

                        infos.curMeshSeedQueue.push(seed2);
                    }
                    neighbour->mark = SEEDED;
                }
            }
            queue.pop();
        }
    }

    void MyAlgorithm::floodComplexGroup(GroupParseInfo& infos, SeedInfo& seed)
    {
        Queue<SeedInfo> q;
        q.emplace(seed);

        while (!q.empty())
        {
            if (q.front().seedFacet->mark != VISITED)
            {
                FH curFace = q.front().seedFacet;
                curFace->mark = VISITED;

                //ItstGraph* ig = new ItstGraph(q.front());
            //    ig->createGraph(q.front());
            //    ig->propInd();

            //    loops.clear();
            //    fh->data->iTri->looplets->getAllCircles(loops);

            //    for (auto& loop : loops)
            //    {
            //        if (ig->classify(loop))
            //            addFacet(loop);
            //    }

            //    for (size_t i = 0; i < 3; i++)
            //    {
            //        ig->getCornerInds(&s, i);
            //        s.fh = ? ;
            //        s.fh->mark = SEEDED;

            //        if (IsSameGroup(s.fh, fh))
            //            q.push(s);
            //        else
            //        {
            //        }
            //    }
            //    delete ig;
            }
            q.pop();
        }
    }

}