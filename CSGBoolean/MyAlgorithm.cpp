#include <unordered_map>
#include <unordered_set>
#include <list>
#include <queue>
#include <deque>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include "MyAlgorithm.h"
#include "ItstAlg.h"
#include "ItstGraph.h"
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

        itst = new ItstAlg(pMeshList);
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

        for (int i = 0; i < 2; i++)
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

    void MyAlgorithm::figureOutFaceInds(VH p, VH q, VH r, int meshId, IIndicatorVector* pinds)
    {
        auto &context = p->data->proxy->pointer()->ctx;
        FH coincident;
        assert(!CGAL::collinear(p->point(), q->point(), r->point()));

        for (auto &ctx : context)
        {
            if (ctx.meshId == meshId) continue;
            Relation rel = relationOfContext(ctx);

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

    void MyAlgorithm::floodComplexGroup(GroupParseInfo& infos, SeedInfo& s)
    {
        Queue<SeedInfo> q;
        auto &ids = s.seedFacet->data->itstTri->meshIds;

        SeedInfo seed;
        seed.seedFacet = s.seedFacet;;
        seed.seedVertex = s.seedVertex;
        seed.indicators.reset(new SampleIndicatorVector(
            *reinterpret_cast<FullIndicatorVector*>(seed.indicators.get()), ids));
        q.emplace(seed);

        Relation *curRelationTable = new Relation[pMeshList->size()];
        for (size_t i = 0; i < pMeshList->size(); i++)
            curRelationTable[i] = static_cast<Relation>((*s.indicators)[i]);

        for (auto& id : ids)
            curRelationTable[id] = REL_NOT_AVAILABLE;

        TestTree testList;
        CSGTreeNode* curTree = copy2(infos.ttree1.get(), infos.curTreeLeaves);
        auto curRelation = ParsingCSGTree((*pMeshList)[infos.curMeshId], curRelationTable,
            pMeshList->size(), curTree, infos.curTreeLeaves, testList);
        assert(!testList.size() || !testList.begin()->testTree->Parent);

        while (!q.empty())
        {
            if (q.front().seedFacet->mark != VISITED)
            {
                auto &curSeed = q.front();
                FH curface = seed.seedFacet;
                curface->mark = VISITED;

                ItstGraph* ig = new ItstGraph(curface, itst, infos.curMeshId);
                assert(ig->get_bValid());

                ig->floodFilling(seed.seedVertex,
                    *reinterpret_cast<SampleIndicatorVector*>(seed.indicators.get()), ids);

                std::deque<ItstGraph::Loop> loops;
                ig->getAllLoops(loops);

                for (auto& loop : loops)
                {
                    if (needAdd(curface, loop, testList))
                        addLoop(loop);
                }

                for (int i = 0; i < 3; i++)
                {
                    auto neighbour = curface->edges[i]->opposite()->facet();

                    SeedInfo seed2;
                    seed2.seedFacet = neighbour;
                    seed2.seedVertex = curface->edges[i]->opposite()->vertex();

                    int gId = curface->edges[i]->opposite()->vertex()->data->proxy->pointer()->idx;
                    SampleIndicatorVector* inds = reinterpret_cast<SampleIndicatorVector*>(ig->get_nodes()[gId].indicator);

                    if (neighbour->mark == SEEDED || neighbour->mark == VISITED) // 跟ray-tracing有关
                        continue;

                    if (isSameGroup(neighbour, curface))
                    {
                        auto sample = new SampleIndicatorVector;
                        *sample = *inds;
                        seed2.indicators.reset(sample);
                        q.push(seed2);
                    }
                    else
                    {
                        auto full = new FullIndicatorVector(*reinterpret_cast<FullIndicatorVector*>(s.indicators.get()));
                        full->fillInSample(inds, ids);
                        seed2.indicators.reset(full);

                        infos.curMeshSeedQueue.push(seed2);
                    }
                    neighbour->mark = SEEDED;
                }

                SAFE_DELETE(ig);
            }
            q.pop();
        }
    }

    bool MyAlgorithm::needAdd(FH fh, ItstGraph::Loop& loop, TestTree& testList)
    {
        // decide a relation
        SampleIndicatorVector sample(fh->data->itstTri->meshIds);
        size_t n = loop.size();
        int checkPoint = -1;

        for (int id : fh->data->itstTri->meshIds)
        {
            Indicator ind = REL_NOT_AVAILABLE;
            int i = 0;
            while (loop[i]->vproxy.pointer()->hasContext(id) && i < n)
                ind = loop[i++]->indicator->at(id);

            if (ind != REL_ON_BOUNDARY)
                sample[id] = ind;
            else
            {
                if (checkPoint == -1)
                {
                    for (size_t i = 0; i < n; i++)
                        if (same_orientation(loop[i]->vproxy.pointer()->pos,
                            loop[i + 1]->vproxy.pointer()->pos,
                            loop[i + 2]->vproxy.pointer()->pos,
                            fh->data->sp.orthogonal_vector()))
                        {
                            checkPoint = i;
                            break;
                        }
                    sample[id] = determineRelationOfFacet(*loop[checkPoint]->vproxy.pointer()->findInContext(id),
                        loop[i + 1]->vproxy.pointer()->pos,
                        loop[i + 1]->vproxy.pointer()->pos);
                }
            }
        }

        // decide retain or drop
        CSGTreeNode* curNode;
        Relation curRelation(REL_UNKNOWN), outRelation(REL_UNKNOWN);
        unsigned testId;
        bool pass = true;
        for (auto &test : testList)
        {
            curNode = GetFirstNode(test.testTree);
            while (curNode)
            {
                testId = curNode->pMesh->Id;
                curRelation = static_cast<Relation>(sample.at(testId));
                curNode = GetNextNode(curNode, curRelation, outRelation);
            }
            if (!(test.targetRelation & outRelation))
            {
                pass = false;
                break;
            }
        }
        return pass;
    }


}