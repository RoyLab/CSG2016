#include "precompile.h"
#include <vector>
#include <memory>
#include <queue>
#include <map>
#include <set>
#include <xlogger.h>

#include "RegularMesh.h"
#include "Octree.h"
#include "intersection.h"
#include "csg.h"
#include "BSP.h"

namespace Boolean
{
    typedef uint32_t MeshId;

    static RegularMesh* g_debug_mesh;

    struct SSeed
    {
        EdgeIndex edgeId;
        IPolygon* pFace = nullptr;
        AutoPtr<IIndicatorVector> eIndicators;

        SSeed() {}
        SSeed(const SSeed& other) { *this = other; }
        SSeed& operator=(const SSeed& other);
    };

    void calcEdgeIndicator(VertexIndex seedVertexId, EdgeIndex seedEdgeId, 
        FullIndicatorVector& vInds, FullIndicatorVector& eInds)
    {
        for (size_t i = 0; i < vInds.getNumber(); i++)
            eInds[i] = vInds[i];

        XLOG_ERROR << UNIMPLEMENTED_DECLARATION;
    }

    void calcNormalSeed(VertexIndex seedId, SSeed& seed, size_t nMesh, Octree* pOctree)
    {
        FullIndicatorVector vInds(nMesh);
        MyVertex& seedV = xvertex(seedId);
        cyPointT seedPts = seedV.vertex_rep();
        for (int i = 0; i < nMesh; i++)
        {
            vInds[i] = REL_UNKNOWN;
        }

        // 找到那些on boundary的mesh，通过遍历所有的边上的面，查看它们的meshId
        for (auto &edgeId : seedV.edges_local())
        {
            MyEdge& eRef = xedge(edgeId);
            MyVertex& theOther = eRef.theOtherVertex(seedId);
            if (theOther.isPlaneRep()) continue; // 过滤掉非原始边

            auto fItr = MyEdge::ConstFaceIterator(eRef);
            for (; fItr; ++fItr)
            {
                if (!fItr.face()) break;
                if (fItr.face()->getType() == IPolygon::SUBPOLYGON) continue; // 过滤掉非原始面
                vInds[fItr.face()->meshId()] = REL_ON_BOUNDARY;
            }
        }

        for (int i = 0; i < nMesh; i++)
        {
            if (vInds[i] == REL_UNKNOWN)
            {
                vInds[i] = PolyhedralInclusionTest(
                    seedPts, 
                    pOctree, 
                    xmeshlist(), 
                    i, 
                    xmeshlist()[i]->inverse()
                );
            }
        }

        calcEdgeIndicator(seedId, seed.edgeId, vInds,
            *(FullIndicatorVector*)seed.eIndicators.get());
    }

    void calcFirstSeed(VertexIndex seedId, SSeed& seed, size_t nMesh)
    {
        FullIndicatorVector vInds(nMesh);
        MyVertex& seedV = xvertex(seedId);
        for (int i = 0; i < nMesh; i++)
            vInds[i] = REL_OUTSIDE;

        // 找到那些on boundary的mesh，通过遍历所有的边上的面，查看它们的meshId
        for (auto &edgeId : seedV.edges_local())
        {
            MyEdge& eRef = xedge(edgeId);
            MyVertex& theOther = eRef.theOtherVertex(seedId);
            if (theOther.isPlaneRep()) continue; // 过滤掉非原始边

            auto fItr = MyEdge::ConstFaceIterator(eRef);
            for (; fItr; ++fItr)
            {
                if (!fItr.face()) break;
                if (fItr.face()->getType() == IPolygon::SUBPOLYGON) continue; // 过滤掉非原始面
                vInds[fItr.face()->meshId()] = REL_ON_BOUNDARY;
            }
        }

        // 找到一个包含有效面，且周围面的不共面的边，赋值edgeId
        // 这里会有隐含的假设：这个极点必须是原始点（应该是对的）
        bool flag = false;
        for (auto &edgeId : seedV.edges_local())
        {
            MyEdge& eRef = xedge(edgeId);
            XPlane base;
            auto fItr = MyEdge::FaceIterator(eRef);
            for (; fItr; ++fItr)
            {
                if (!fItr.face()->isValid()) continue;

                if (fItr.face()->getType() == IPolygon::TRIANGLE)
                {
                    Triangle* pTri = (Triangle*)fItr.face();
                    try
                    {
                        pTri->calcSupportingPlane();
                    }
                    catch (...)
                    {
                        XLOG_ERROR << "degenerate triangle";
                        pTri->mark = VISITED;
                        continue;
                    }
                }
                
                if (!base.is_valid())
                {
                    base = fItr.face()->supportingPlane();
                }
                else
                {
                    if (!base.parallel(fItr.face()->supportingPlane()))
                    {
                        seed.edgeId = edgeId;
                        flag = true;
                        break;
                    }
                }
            }
            if (flag) break;
        }
        assert(vertex_id_equals_simple(xedge(seed.edgeId).ends[0], seedId)
            || vertex_id_equals_simple(xedge(seed.edgeId).ends[1] ,seedId));

        if (!flag)
        {
            throw "";
        }

        // 赋值pFace
        MyEdge& eRef = xedge(seed.edgeId);
        auto fItr = MyEdge::FaceIterator(eRef);
        flag = false;
        for (; fItr; ++fItr)
        {
            if (fItr.face()->isValid())
            {
                flag = true;
                seed.pFace = fItr.face();
                break;
            }
        }
        assert(flag);

        calcEdgeIndicator(seedId, seed.edgeId, vInds, 
            *(FullIndicatorVector*)seed.eIndicators.get());
    }

    Relation vRelation2fRelation(Oriented_side rel, XPlane& testPlane, XPlane& refPlane)
    {
        if (rel == ON_POSITIVE_SIDE)
        return REL_OUTSIDE;
        else if (rel == ON_NEGATIVE_SIDE)
            return REL_INSIDE;
        else
        {
            const Real* dataA = testPlane.get_data();
            const Real* dataB = refPlane.get_data();

            if (dataA[0] * dataB[0] > 0 || dataA[1] * dataB[1] > 0
                || dataA[2] * dataB[2] > 0)
                return REL_SAME;
            else
                return REL_OPPOSITE;
        }
    }

    void calcFaceIndicator(SSeed& seed, std::vector<Relation>& relTab)
    {
        MyEdge& edge = xedge(seed.edgeId);
        IPolygon* polygon = seed.pFace;

        // copy the relation
        const size_t nMesh = ((FullIndicatorVector*)seed.eIndicators.get())->getNumber();
        for (size_t i = 0; i < nMesh; i++)
        {
            relTab[i] = (Relation)seed.eIndicators.get()->at(i);
//#ifdef XR_DEBUG
//            bool flag = true;
//            if (i != polygon->meshId() && relTab[i] == REL_ON_BOUNDARY)
//            {
//                flag = false;
//                for (auto& neigh : *edge.neighbor)
//                {
//                    if (neigh.first == i)
//                    {
//                        flag = true;
//                        break;
//                    }
//                }
//            }
//#endif
        }

        relTab[polygon->meshId()] = REL_SAME;

        if (!edge.neighbor || edge.neighbor->empty()) return;

        MyVertex& repVertex = xvertex(seed.pFace->get_rep_vertex(seed.edgeId));

        // correct the relation
        for (auto &neigh: *edge.neighbor)
        {
            Oriented_side side;

            if (neigh.first == polygon->meshId())
                continue;

            if (neigh.second.type == NeighborInfo::Edge)
            {
                std::vector<Triangle*> faces;
                auto fItr = MyEdge::FaceIterator(xedge(neigh.second.neighborEdgeId), true);
                for (; fItr && fItr.face(); fItr.incrementToTriangle())
                {
                    faces.push_back(fItr.as_triangle());
                }

                BSPTree bsp; XPlane bspPlane;
                bsp.buildNoCross(faces);
                side = bsp.classify(repVertex, &bspPlane);

                relTab[neigh.first] = vRelation2fRelation(side,
                    bspPlane, polygon->supportingPlane());
            }
            else
            {
                assert(neigh.second.type == NeighborInfo::Face);
                Oriented_side side = orientation(neigh.second.pTrangle->supportingPlane(), repVertex);

                relTab[neigh.first] = vRelation2fRelation(side,
                    neigh.second.pTrangle->supportingPlane(), polygon->supportingPlane());
            }
        }
    }

    void doClassification(Octree* pOctree, CSGTree<RegularMesh>* pCSG,
        std::vector<RegularMesh*>& meshList, RegularMesh* result,
        VertexIndex seedId, RegularMesh* debug_mesh)
    {
        g_debug_mesh = debug_mesh;
        GlobalData* pMem = GlobalData::getObject();

        uint32_t nMesh = meshList.size();
        CSGTreeOld* tree = pCSG->auxiliary();
        CSGTreeNode** curTreeLeaves = new CSGTreeNode*[nMesh];

        SSeed tmpSeed; tmpSeed.eIndicators.reset(new FullIndicatorVector(nMesh));
        MeshId curMeshId = 0;
        std::queue<IPolygon*> faceQueue;
        std::queue<SSeed> intraQueue, interQueue;
        bool added = false, inverse = false;
        std::pair<uint32_t, uint32_t> inversePair;
        Relation relation = REL_NOT_AVAILABLE;
        TestTree dummyForest;
        std::vector<Relation> relTab(nMesh);
        std::vector<EdgeIndex> edges;

        int seed_count = 0;
        for (int imesh = 0; imesh < xmeshlist().size(); ++imesh)
        {
            RegularMesh* mesh = xmeshlist()[imesh];
            for (int iface = 0; iface < mesh->faces().size(); ++iface)
            {
                IPolygon* cur_face = mesh->faces()[iface];
                if (!cur_face->isValid() ||
                    cur_face->mark == VISITED)
                {
                    continue;
                }

                if (seed_count == 0)
                {
                    MyVertex& seedV = xvertex(seedId);
                    tmpSeed.edgeId = *seedV.edges_local().begin();
                    calcFirstSeed(seedId, tmpSeed, nMesh);
                    ++seed_count;
                }
                else
                {
                    bool find_orig = false;
                    int i;
                    for (i = 0; i < cur_face->outer_degree(); ++i)
                    {
                        if (pMem->is_original_vertex(cur_face->outer_vertex_id(i)))
                        {
                            find_orig = true;
                            break;
                        }
                    }

                    if (!find_orig)
                    {
                        continue;
                    }

                    std::vector<EdgeIndex> edges;
                    cur_face->getAllEdges(edges);

                    tmpSeed.edgeId = edges[0];
                    tmpSeed.pFace = cur_face;
                    calcNormalSeed(cur_face->outer_vertex_id(i), tmpSeed, nMesh, pOctree);

                    XLOG_DEBUG << "new seed " << seed_count;
                }

                interQueue.push(tmpSeed);
                tmpSeed.pFace->mark = SEEDED0;
                assert(tmpSeed.pFace->isValid());
                while (!interQueue.empty())
                {
                    SSeed curSeed = interQueue.front();
                    interQueue.pop();

                    if (curSeed.pFace->mark == VISITED)
                        continue;

                    assert(intraQueue.empty());
                    curSeed.pFace->mark = SEEDED1;
                    intraQueue.push(curSeed);

                    curMeshId = curSeed.pFace->meshId();
                    inverse = meshList[curMeshId]->inverse();
                    if (inverse)
                        inversePair.first = result->faces().size();

                    //bool seedFlag = true;
                    while (!intraQueue.empty())
                    {
                        curSeed = intraQueue.front();
                        intraQueue.pop();

                        if (curSeed.pFace->mark == VISITED)
                            continue;

                        CSGTreeNode* tree0 = copy2(tree->pRoot, curTreeLeaves);

                        try
                        {
                            calcFaceIndicator(curSeed, relTab);
                        }
                        catch (int e)
                        {
                            if (e == 1)
                            {
                                curSeed.pFace->mark = VISITED;
                                continue;
                            }
                        }
                        //if (seedFlag) seedFlag = false; // only for debug use

                        relation = ParsingCSGTree(meshList[curMeshId], relTab.data(),
                            nMesh, tree0, curTreeLeaves, dummyForest);
                        assert(relation != REL_NOT_AVAILABLE || relation != REL_UNKNOWN);

                        added = (relation == REL_SAME);
                        faceQueue.push(curSeed.pFace);
                        while (!faceQueue.empty())
                        {
                            IPolygon* curFace = faceQueue.front();
                            faceQueue.pop();
                            if (curFace->mark == VISITED)
                                continue;

                            curFace->mark = VISITED;
                            if (added)
                                result->faces().push_back(curFace);

                            edges.clear();
                            curFace->getAllEdges(edges);
                            //assert(edges.size() == curFace->degree());

                            // flood filling, bfs
                            for (int i = 0; i < edges.size(); i++)
                            {
                                MyEdge& curEdge = xedge(edges[i]);
                                MyEdge::FaceIterator fItr(curEdge);
                                if (curEdge.neighbor)
                                {
                                    // 但因为共面的存在，有些neigh可能记录的相邻，但不在相邻面当中，这没关系！
                                    //assert(curEdge.faceCount() >= 4); //  ，那么一定会有超过4个polygon在周围
                                    for (; fItr; ++fItr)
                                    {
                                        if (!fItr.face())
                                        {
                                            XLOG_ERROR << "Edge with less than two neighboring faces.";
                                            continue;
                                        }

                                        if (!fItr.face()->isValid()) continue;

                                        tmpSeed.edgeId = edges[i];
                                        tmpSeed.pFace = fItr.face();

                                        for (int i = 0; i < nMesh; i++)
                                        {
                                            Indicator tmpInd = REL_NOT_AVAILABLE;
                                            switch (relTab[i])
                                            {
                                            case REL_SAME:
                                            case REL_OPPOSITE:
                                                tmpInd = REL_ON_BOUNDARY;
                                                break;
                                            default:
                                                tmpInd = relTab[i];
                                                break;
                                            }
                                            tmpSeed.eIndicators->at(i) = tmpInd;
                                        }

                                        for (auto& nInfo : *curEdge.neighbor)
                                            tmpSeed.eIndicators->at(nInfo.first) = REL_ON_BOUNDARY;

                                        if (fItr.face()->meshId() == curMeshId)
                                        {
                                            if (fItr.face()->mark < SEEDED1)
                                            {
                                                tmpSeed.pFace->mark = SEEDED1;
                                                intraQueue.push(tmpSeed);
                                            }
                                        }
                                        else
                                        {
                                            if (fItr.face()->mark < SEEDED0)
                                            {
                                                tmpSeed.pFace->mark = SEEDED0;
                                                interQueue.push(tmpSeed);
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    for (; fItr; ++fItr)
                                    {
                                        if (!fItr.face())
                                        {
                                            XLOG_ERROR << "Edge with less than two neighboring faces.";
                                            continue;
                                        }
                                        if (!fItr.face()->isValid() || fItr.face()->meshId() != curMeshId) continue;
                                        if (fItr.face()->mark < SEEDED2)
                                        {
                                            faceQueue.push(fItr.face());
                                            fItr.face()->mark = SEEDED2;
                                        }
                                    }
                                }
                            }
                        }
                    } // loop inseed mesh

                    if (inverse)
                    {
                        inversePair.second = result->faces().size();
                        result->inverseMap.push_back(inversePair);
                    }
                } // loop seed
            } // face loop
        } // mesh loop

        SAFE_DELETE_ARRAY(curTreeLeaves);
    }

    SSeed & SSeed::operator=(const SSeed & other)
    {
        edgeId = other.edgeId;
        pFace = other.pFace;

        if (other.eIndicators->getType() == IIndicatorVector::FULL)
            eIndicators.reset(new FullIndicatorVector(
                *reinterpret_cast<FullIndicatorVector*>(other.eIndicators.get())));
        else eIndicators.reset(new SampleIndicatorVector(
                *reinterpret_cast<SampleIndicatorVector*>(other.eIndicators.get())));

        return *this;
    }

}