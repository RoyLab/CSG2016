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

    struct SSeed
    {
        MyEdge::Index edgeId;
        IPolygon* pFace = nullptr;
        AutoPtr<IIndicatorVector> eIndicators;

        SSeed() {}
        SSeed(const SSeed& other) { *this = other; }
        SSeed& operator=(const SSeed& other);
    };

    void calcEdgeIndicator(MyVertex::Index seedVertexId, MyEdge::Index seedEdgeId, 
        FullIndicatorVector& vInds, FullIndicatorVector& eInds)
    {
        for (size_t i = 0; i < vInds.getNumber(); i++)
            eInds[i] = vInds[i];

        XLOG_ERROR << UNIMPLEMENTED_DECLARATION;
    }

    void calcEdgeIndicatorByExtremity(MyVertex::Index seedId, SSeed& seed, FullIndicatorVector& inds, size_t nMesh)
    {
        FullIndicatorVector vInds(nMesh);
        MyVertex& seedV = xvertex(seedId);
        for (int i = 0; i < nMesh; i++)
            vInds[i] = REL_OUTSIDE;

        // 找到那些on boundary的mesh，通过遍历所有的边上的面，查看它们的meshId
        for (auto &edgeId : seedV.edges)
        {
            MyEdge& eRef = xedge(edgeId);
            MyVertex& theOther = eRef.theOtherVertex(seedId);
            if (theOther.isPlaneRep()) continue; // 过滤掉非原始边

            auto fItr = MyEdge::ConstFaceIterator(eRef);
            for (; fItr; ++fItr)
            {
                if (fItr.face()->getType() == IPolygon::SUBPOLYGON) continue; // 过滤掉非原始面
                vInds[fItr.face()->meshId()] = REL_ON_BOUNDARY;
            }
        }

        // 找到一个包含有效面，且周围面的不共面的边，赋值edgeId
        // 这里会有隐含的假设：这个极点必须是原始点（应该是对的）
        bool flag = false;
        for (auto &edgeId : seedV.edges)
        {
            MyEdge& eRef = xedge(edgeId);
            MyVertex& theOther = eRef.theOtherVertex(seedId);

            auto fItr = MyEdge::FaceIterator(eRef);
            if (fItr.face()->getType() == IPolygon::TRIANGLE)
            {
                reinterpret_cast<Triangle*>(fItr.face())->calcSupportingPlane();
            }
            XPlane basePlane = fItr.face()->supportingPlane();
            flag = false; // 是否有不共面的相邻面
            for (; fItr; ++fItr)
            {
                if (!fItr.face()->isValid()) continue;

                if (fItr.face()->getType() == IPolygon::TRIANGLE)
                {
                    Triangle* pTri = (Triangle*)fItr.face();
                    pTri->calcSupportingPlane();
                }
                
                if (!basePlane.coplanar(fItr.face()->supportingPlane()))
                {
                    seed.edgeId = edgeId;
                    flag = true;
                    break;
                }
            }
            if (flag) break;
        }
        assert(flag);
        assert(xedge(seed.edgeId).ends[0] == seedId || xedge(seed.edgeId).ends[1] == seedId);

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

        calcEdgeIndicator(seedId, seed.edgeId, vInds, inds);
    }

    Relation vRelation2fRelation(Oriented_side rel, XPlane& testPlane, XPlane& refPlane)
    {
        if (rel == ON_POSITIVE_SIDE)
        return REL_OUTSIDE;
        else if (rel == ON_NEGATIVE_SIDE)
            return REL_INSIDE;
        else
        {
            const Real* dataA = testPlane.data();
            const Real* dataB = refPlane.data();

            if (dataA[0] * dataB[0] > 0 || dataA[1] * dataB[1] > 0
                || dataA[2] * dataB[2] > 0)
                return REL_SAME;
            else
                return REL_OPPOSITE;
        }
    }

    MyVertex::Index GetRepVertex_SPOLY(MyEdge::Index edgeId, SubPolygon* polygon)
    {
        assert(orientation(polygon->supportingPlane(),(polygon->vertex(0))) == ON_ORIENTED_BOUNDARY);
        assert(orientation(polygon->supportingPlane(),(polygon->vertex(1))) == ON_ORIENTED_BOUNDARY);
        assert(orientation(polygon->supportingPlane(),(polygon->vertex(2))) == ON_ORIENTED_BOUNDARY);

        MyEdge& edge = xedge(edgeId);

        // find init point
        int edgeIndexInFace = -1;
        for (int i = 0; i < polygon->degree(); i++)
        {
            if (polygon->edgeId(i) == edgeId)
            {
                edgeIndexInFace = i;
                break;
            }
        }
        assert(edgeIndexInFace != -1);

        MyVertex::Index vIdInPlane;
        XPlane boundPlane, tmpPlane;
        for (int i = 2; i < polygon->degree(); i++)
        {
            vIdInPlane = polygon->vertexId((edgeIndexInFace + i) % polygon->degree());

            // find a bounding plane
            bool flag = false;
            for (auto &neigh : *edge.neighbor)
            {
                if (neigh.type == NeighborInfo::Edge)
                {
                    for (auto fItr = MyEdge::FaceIterator(xedge(neigh.neighborEdgeId), true);
                        fItr; fItr.incrementToTriangle())
                    {
                        assert(fItr.face()->getType() == IPolygon::TRIANGLE);
                        tmpPlane = ((Triangle*)fItr.face())->supportingPlane();
                        if (orientation(tmpPlane, xvertex(vIdInPlane)) != ON_ORIENTED_BOUNDARY)
                        {
                            boundPlane = tmpPlane;
                            flag = true;
                            break;
                        }
                    }
                }
                else
                {
                    assert(neigh.type == NeighborInfo::Face);
                    XPlane tmpPlane = neigh.pTrangle->supportingPlane();
                    if (orientation(tmpPlane, xvertex(vIdInPlane)) != ON_ORIENTED_BOUNDARY)
                    {
                        boundPlane = tmpPlane;
                        break;
                    }
                }
                if (flag) break;
            }
            if (boundPlane.isValid()) break;
        }
        assert(boundPlane.isValid());

        // correct the direction of bounding plane
        XLine edgeLine(polygon->supportingPlane(), boundPlane);
        assert(!polygon->supportingPlane().idEquals(boundPlane));
        int tmpSide = linearOrder(edgeLine, xvertex(polygon->vertexId((edgeIndexInFace + 1) % polygon->degree())),
            xvertex(polygon->vertexId(edgeIndexInFace)));

        assert(tmpSide != 0);
        if (tmpSide < 0)
            boundPlane.inverse();

        // pick a correct rep vertex
        MyVertex::Index repVertexId = INVALID_UINT32;
        for (size_t i = 0; i < polygon->degree(); i++)
        {
            if (orientation(boundPlane, xvertex(polygon->vertexId(i))) == ON_POSITIVE_SIDE)
            {
                repVertexId = polygon->vertexId(i);
                break;
            }
        }
        XR_assert(repVertexId != INVALID_UINT32);
        return repVertexId;
    }

    MyVertex::Index GetRepVertex_TRI(MyEdge::Index edgeId, Triangle* polygon)
    {
        MyEdge& edge = xedge(edgeId);

        int edgeIndexInFace = -1;
        for (size_t i = 0; i < polygon->degree(); i++)
        {
            if (polygon->edgeId(i) == edgeId)
            {
                edgeIndexInFace = i;
                break;
            }
        }
        assert(edgeIndexInFace != -1);
        return polygon->vertexId(edgeIndexInFace);
    }

    void calcFaceIndicator(SSeed& seed, std::vector<Relation>& relTab, bool hasNeighbor)
    {
        MyEdge& edge = xedge(seed.edgeId);
        IPolygon* polygon = seed.pFace;

        // copy the relation
        const size_t nMesh = ((FullIndicatorVector*)seed.eIndicators.get())->getNumber();
        for (size_t i = 0; i < nMesh; i++)
        {
            relTab[i] = (Relation)seed.eIndicators.get()->at(i);
#ifdef XR_DEBUG
            bool flag = true;
            if (i != polygon->meshId() && relTab[i] == REL_ON_BOUNDARY)
            {
                flag = false;
                for (auto& neigh : *edge.neighbor)
                {
                    if (neigh.neighborMeshId == i)
                    {
                        flag = true;
                        break;
                    }
                }
            }
            assert(flag);
#endif
        }

        relTab[polygon->meshId()] = REL_SAME;

        assert(!hasNeighbor || edge.neighbor->size());
        if (!edge.neighbor || edge.neighbor->empty()) return;

        // remove repetitive neighborInfo, must before <find repVertex>
        // because subpolygon use neighborInfo to find a splitPlane
        if (!edge.noOverlapNeighbor)
        {
            edge.noOverlapNeighbor = true;
            std::set<MeshId> meshSets;
            std::vector<NeighborInfo> newNeighbor;
            for (auto &neigh : *edge.neighbor)
            {
                if (meshSets.find(neigh.neighborMeshId) == meshSets.end())
                {
                    newNeighbor.push_back(neigh);
                    meshSets.insert(neigh.neighborMeshId);
                }
            }
            edge.neighbor->swap(newNeighbor);
        }

        // find the represented vertex
        MyVertex::Index repVertexId;
        if (seed.pFace->getType() == IPolygon::TRIANGLE)
            repVertexId = GetRepVertex_TRI(seed.edgeId, (Triangle*)seed.pFace);
        else
            repVertexId = GetRepVertex_SPOLY(seed.edgeId, (SubPolygon*)seed.pFace);
        MyVertex& repVertex = xvertex(repVertexId);

        // correct the relation
        for (auto &neigh: *edge.neighbor)
        {
            Oriented_side side;

            if (neigh.neighborMeshId == polygon->meshId())
                continue;

            if (neigh.type == NeighborInfo::Edge)
            {
                std::vector<IPolygon*> faces;
                auto fItr = MyEdge::FaceIterator(xedge(neigh.neighborEdgeId), true);
                for (; fItr; fItr.incrementToTriangle())
                {
                    assert(fItr.face()->getType() == IPolygon::TRIANGLE);
                    faces.push_back(fItr.face());
                }

                BSPTree bsp; XPlane bspPlane;
                bsp.buildNoCross(faces);
                side = bsp.classify(repVertex, &bspPlane);

                relTab[neigh.neighborMeshId] = vRelation2fRelation(side,
                    bspPlane, polygon->supportingPlane());
            }
            else
            {
                assert(neigh.type == NeighborInfo::Face);
                Oriented_side side = orientation(neigh.pTrangle->supportingPlane(), repVertex);

                relTab[neigh.neighborMeshId] = vRelation2fRelation(side,
                    neigh.pTrangle->supportingPlane(), polygon->supportingPlane());
            }
        }
    }

    void doClassification(Octree* pOctree, CSGTree<RegularMesh>* pCSG,
        std::vector<RegularMesh*>& meshList, RegularMesh* result,
        MyVertex::Index seedId)
    {
        uint32_t nMesh = meshList.size();
        CSGTreeOld* tree = pCSG->auxiliary();
        CSGTreeNode** curTreeLeaves = new CSGTreeNode*[nMesh];

        SSeed tmpSeed;
        MyVertex& seedV = xvertex(seedId);
        tmpSeed.edgeId = *seedV.edges.begin();
        MyEdge& seedE = xedge(tmpSeed.edgeId);

        tmpSeed.eIndicators.reset(new FullIndicatorVector(nMesh));
        calcEdgeIndicatorByExtremity(seedId, tmpSeed, 
            *reinterpret_cast<FullIndicatorVector*>(tmpSeed.eIndicators.get()), nMesh);

        MeshId curMeshId;
        std::queue<IPolygon*> faceQueue;
        std::queue<SSeed> intraQueue, interQueue;
        bool added, inverse;
        std::pair<uint32_t, uint32_t> inversePair;
        Relation relation;
        TestTree dummyForest;
        std::vector<Relation> relTab(nMesh);
        std::vector<MyEdge::Index> edges;

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

            bool seedFlag = true;
            while (!intraQueue.empty())
            {
                curSeed = intraQueue.front();
                intraQueue.pop();

                if (curSeed.pFace->mark == VISITED)
                    continue;

                CSGTreeNode* tree0 = copy2(tree->pRoot, curTreeLeaves);

                calcFaceIndicator(curSeed, relTab, seedFlag?false:true);
                if (seedFlag) seedFlag = false; // only for debug use

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
                    curFace->getEdges(edges);
                    assert(edges.size() == curFace->degree());

                    // flood filling, bfs
                    for (int i = 0; i < curFace->degree(); i++)
                    {
                        MyEdge& curEdge = xedge(edges[i]);
                        MyEdge::FaceIterator fItr(curEdge);
                        if (curEdge.neighbor)
                        {
                            // 但因为共面的存在，有些neigh可能记录的相邻，但不在相邻面当中，这没关系！
                            //assert(curEdge.faceCount() >= 4); // 如果这是一个相交而成的边，那么一定会有超过4个polygon在周围
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

                                for (NeighborInfo& nInfo : *curEdge.neighbor)
                                    tmpSeed.eIndicators->at(nInfo.neighborMeshId) = REL_ON_BOUNDARY;

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
                //break;
            }
            if (inverse)
            {
                inversePair.second = result->faces().size();
                result->inverseMap.push_back(inversePair);
            }
            //break;
        }
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