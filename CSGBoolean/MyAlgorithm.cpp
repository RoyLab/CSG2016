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
    inline bool getSharedEdge(FH f0, FH f1, EH& eh)
    {
        for (size_t i = 0; i < 3; i++)
        {
            if (f0->edges[i]->opposite()->facet() == f1)
            {
                eh = f0->edges[i];
                return true;
            }
        }
        return false;
    }

    void MyAlgorithm::solve(const std::string& expr, std::vector<MyMesh*>& meshes)
    {
        pMeshList = &meshes;
        for (auto mesh : *pMeshList)
            mesh->calcBbox();

        CGAL::Bbox_3 aabb = pMeshList->at(0)->get_bbox();
        for (auto mesh : *pMeshList)
            aabb += mesh->get_bbox();

        for (auto mesh : *pMeshList)
        {
            mesh->normalize(aabb);
            mesh->filter();
            mesh->init();
        }

        CSGTree<MyMesh>* pCsg = new CSGTree<MyMesh>;
        pCsg->createCSGTreeFromExpr(expr, pMeshList->data(), pMeshList->size());
        pCsg->makePositiveAndLeftHeavy();

        pOctree = new Octree;
        std::vector<Octree::Node*> intersectLeaves;
        pOctree->build(*pMeshList, &intersectLeaves);

        itst = new ItstAlg(pMeshList);
        itst->doIntersection(intersectLeaves);

        std::cout.precision(18);
        for (auto pt : itst->vEnt)
            std::cout << "point: " << pt->pos.getCoord() << std::endl;

        floodColoring(pCsg, itst);

        csgResult->denormalize(aabb);

        SAFE_DELETE(itst);
        SAFE_DELETE(pCsg);
        SAFE_DELETE(pOctree);
    }

    IIndicatorVector* MyAlgorithm::computeFullIndicator(VH vh, size_t meshId)
    {
        IIndicatorVector *pind = new FullIndicatorVector(pMeshList->size());
        IIndicatorVector &ind = *pind;
        ind[meshId] = REL_ON_BOUNDARY;

        if (vh->data && vh->data->proxy)
        {
            const auto& ctx = (*vh->data->proxy).pointer()->ctx;
            for (auto &c : ctx)
                ind[c.meshId] = REL_ON_BOUNDARY;
        }

        for (size_t meshId = 0; meshId < pMeshList->size(); meshId++)
        {
            if (ind[meshId] == REL_ON_BOUNDARY)
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
                std::vector<int> itstPrims;
                itstAlg->get_adjGraph()->getIntersectPrimitives(curSeed.meshId, itstPrims);

                CSGTreeNode* tree0 = copy2(tree->pRoot, infos.curTreeLeaves);
                for (size_t id = 0; id < nMesh; id++)
                    infos.curTreeLeaves[id]->relation = static_cast<Relation>(curSeed.indicators->at(id));

                for (size_t id : itstPrims)
                    infos.curTreeLeaves[id]->relation = REL_UNKNOWN;

                Relation meshRel = CompressCSGNodeIteration(tree0);
                infos.ttree1.reset(tree0);

                if (meshRel == REL_NOT_AVAILABLE)
                {
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
                        }
                        infos.curMeshSeedQueue.pop();
                    }
                }
                infos.otherMeshSeedQueue.pop();
            }
        }
        csgResult->delegate(*pConstruct);
        SAFE_DELETE(pConstruct);
        SAFE_DELETE_ARRAY(infos.curTreeLeaves);
    }

    void MyAlgorithm::figureOutFaceInds(SeedInfo& s, int meshId, Relation** dirRelation)
    {
        if (s.seedVertex->data && s.seedVertex->data->proxy->pointer()->ctx.size() > 1)
        {
            assert(s.seedVertex->halfedge()->facet() == s.seedFacet);
            EH eh = s.seedVertex->halfedge();

            assert(eh->vertex() == s.seedVertex);
            VH pts[2];
            for (size_t i = 0; i < 2; i++)
            {
                eh = eh->next();
                pts[i] = eh->vertex();
            }

            auto &p = s.seedVertex, &q = pts[0], &r = pts[1];
            assert(!CGAL::collinear(p->point(), q->point(), r->point()));
            auto ventity = p->data->proxy->pointer();

            for (auto &ctx : ventity->ctx)
            {
                Relation relation = REL_UNKNOWN;
                if (ctx.meshId == meshId)
                    relation = REL_SAME;
                else
                    relation = determineRelationOfFacet(ctx, q->data->proxy->pointer()->pos,
                        r->data->proxy->pointer()->pos, s.seedFacet->data->sp.orthogonal_vector());

                (*dirRelation)[ctx.meshId] = relation;
            }
        }

        (*dirRelation)[meshId] = REL_SAME;
        for (size_t i = 0; i < pMeshList->size(); i++)
        {
            assert((*dirRelation)[i] != REL_ON_BOUNDARY);
            if ((*dirRelation)[i] == REL_SAME || (*dirRelation)[i] == REL_OPPOSITE)
                s.indicators->at(i) = REL_ON_BOUNDARY;
            else
                s.indicators->at(i) = (*dirRelation)[i];
        }
    }

    void MyAlgorithm::AddFacet(FH fh)
    {
        for (int i = 0; i < 3; i++)
        {
            if (fh->vertices[i]->idx == -1)
            {
                if (fh->vertices[i]->data && fh->vertices[i]->data->proxy->pointer()->resultId != -1)
                    fh->vertices[i]->idx = fh->vertices[i]->data->proxy->pointer()->resultId;
                else
                {
                    int idx = pConstruct->add_vertex(fh->vertices[i]->point());
                    fh->vertices[i]->idx = idx;
                    if (fh->vertices[i]->data)
                        fh->vertices[i]->data->proxy->pointer()->resultId = idx;
                }
            }
        }

        pConstruct->addFacets(fh->vertices[0]->idx,
            fh->vertices[1]->idx, fh->vertices[2]->idx);
    }

    bool MyAlgorithm::isSameGroup(FH fh0, FH fh1) const
    {
        assert(fh0 != fh1);
        if (fh0->isSimple() ^ fh1->isSimple())
            return false;

        if (fh0->isSimple())
            return true;

        return fh0->data->itstTri->meshIds == fh1->data->itstTri->meshIds;
    }

    void MyAlgorithm::genVertexInds(IIndicatorVector* target, VH vh) const
    {
        auto &t = *target;

        if (!vh->data || vh->data->proxy->pointer()->ctx.size() == 0)
            return;

        for (auto& ctx : vh->data->proxy->pointer()->ctx)
            t[ctx.meshId] = REL_ON_BOUNDARY;
    }

    void MyAlgorithm::floodSimpleGroup(GroupParseInfo& infos, SeedInfo& s)
    {
        // 找到这一群的共同indicators
        Relation *curRelationTable = new Relation[pMeshList->size()];
        for (size_t i = 0; i < pMeshList->size(); i++)
            curRelationTable[i] = static_cast<Relation>((*s.indicators)[i]);

        figureOutFaceInds(s, infos.curMeshId, &curRelationTable);

        Queue<FH> queue;
        queue.push(s.seedFacet);

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
            *reinterpret_cast<FullIndicatorVector*>(s.indicators.get()), ids));
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
                FH curface = curSeed.seedFacet;
                curface->mark = VISITED;

                static int count = 0;
                count++;

                ItstGraph* ig = new ItstGraph(curface, itst, infos.curMeshId);
                assert(ig->get_bValid());

                ig->floodFilling(curSeed.seedVertex,
                    *reinterpret_cast<SampleIndicatorVector*>(curSeed.indicators.get()), ids);

                std::deque<ItstGraph::Loop> loops;
                ig->getAllLoops(loops);

                for (auto& loop : loops)
                {
                    if (needAdd(curface, loop, testList))
                        addLoop(loop);
                }

                for (int i = 0; i < 3; i++)
                {
                    auto halfedge = curface->edges[i]->opposite();
                    auto neighbor = halfedge->facet();

                    SeedInfo seed2;
                    seed2.seedFacet = halfedge->face();
                    seed2.seedVertex = halfedge->vertex();

                    int gId = seed2.seedVertex->data->proxy->pointer()->idx; 
                    int localId = ig->get_maps()[gId];
                    SampleIndicatorVector* inds = reinterpret_cast<SampleIndicatorVector*>(ig->get_nodes()[localId].indicator);

                    if (neighbor->mark == SEEDED || neighbor->mark == VISITED) // 跟ray-tracing有关
                        continue;

                    if (isSameGroup(neighbor, curface))
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
                    neighbor->mark = SEEDED;
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
            for (size_t i = 0; i < n; i++)
            {
                if (!loop[i]->vproxy.pointer()->hasContext(id))
                {
                    ind = loop[i]->indicator->at(id);
                    break;
                }
            }

            if (ind != REL_NOT_AVAILABLE)
            {
                sample[id] = ind;
                if (ind == REL_ON_BOUNDARY)
                {
                    ReportError("Not implemented!");
                }
            }
            else
            {
                if (checkPoint == -1)
                {
                    for (size_t j = 0; j < n; j++)
                        if (same_orientation(loop[j]->vproxy.pointer()->pos,
                            loop[j+1]->vproxy.pointer()->pos,
                            loop[j+2]->vproxy.pointer()->pos,
                            fh->data->sp.orthogonal_vector()))
                        {
                            checkPoint = j;
                            break;
                        }
                    sample[id] = determineRelationOfFacet(*loop[checkPoint]->vproxy.pointer()->findInContext(id),
                        loop[checkPoint + 1]->vproxy.pointer()->pos,
                        loop[checkPoint + 2]->vproxy.pointer()->pos, fh->data->sp.orthogonal_vector());
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
                assert(sample[testId] != REL_NOT_AVAILABLE && sample[testId] != REL_UNKNOWN);
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

    void MyAlgorithm::addLoop(ItstGraph::Loop& loop)
    {
        int *indices = new int[loop.size()];
        for (int i = 0; i < loop.size(); i++)
        {
            auto pV = loop[i]->vproxy.pointer();
            assert(pV);
            if (pV->resultId == -1)
            {
                int idx = pConstruct->add_vertex(pV->pos.getCoord());
                pV->resultId = idx;
            }
            indices[i] = pV->resultId;
        }

        pConstruct->addSurface(indices, loop.size());
    }


    //int baoshu = 1;
    //int zuobiao = 0;
    //int renshu = renshu0;
    //while (renshu >１)
    //{
    //    if (!arr[zuobiao].isOut)
    //    {
    //        if (baoshu == 3)
    //        {
    //            arr[zuobiao].isOut = true;
    //            renshu--;
    //        }

    //        baoshu += 1;
    //        if (baoshu == 4)
    //            baoshu = 1;
    //    }
    //    zuobiao = (zuobiao + 1) % renshu0;
    //}
}