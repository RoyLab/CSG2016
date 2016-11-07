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
    template <class T> struct AutoPtr : std::shared_ptr<T> {};
    const int MARK_BEGIN = 0xff; // 因为mark还用来在第三阶段标志有没有被访问过，所以这里让出256个数字用于这些工作????
    enum Mark { UNVISITED, SEEDED0, SEEDED1, SEEDED2, VISITED };
    typedef uint32_t MeshId;

    int linearOrder(const XLine& l, const MyVertex& a, const MyVertex& b);
    struct SSeed
    {
        MyEdge::Index edgeId;
        IPolygon* pFace = nullptr;
        AutoPtr<IIndicatorVector> eIndicators;

        SSeed() {}
        SSeed(const SSeed& other) { *this = other; }
        SSeed& operator=(const SSeed& other);
    };

    Oriented_side orientation(const XPlane& p, const MyVertex& v)
    {
        if (v.isPlaneRep())
            return p.orientation(v.ppoint());
        else return p.orientation(v.point());
    }

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

        for (auto &edgeId : seedV.edges)
        {
            MyEdge& eRef = xedge(edgeId);
            MyVertex& theOther = eRef.theOtherVertex(seedId);
            if (theOther.isPlaneRep()) continue;

            auto fItr = MyEdge::ConstFaceIterator(eRef);
            for (; fItr; ++fItr)
            {
                if (fItr.face()->getType() == IPolygon::SUBPOLYGON) continue;
                vInds[fItr.face()->meshId()] = REL_ON_BOUNDARY;
            }
        }

        bool flag = false;
        for (auto &edgeId : seedV.edges)
        {
            MyEdge& eRef = xedge(edgeId);
            MyVertex& theOther = eRef.theOtherVertex(seedId);

            auto fItr = MyEdge::FaceIterator(eRef);
            for (; fItr; ++fItr)
            {
                if (fItr.face()->isValid())
                {
                    seed.edgeId = edgeId;
                    seed.pFace = fItr.face();
                    flag = true;
                    break;
                }
            }
            if (flag) break;
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
                || dataA[1] * dataB[1] > 0)
                return REL_SAME;
            else
                return REL_OPPOSITE;
        }
    }

    void calcFaceIndicator(SSeed& seed, std::vector<Relation>& relTab)
    {
        MyEdge& edge = xedge(seed.edgeId);
        IPolygon* polygon = seed.pFace;
        assert(edge.neighbor);

        int edgeIndexInFace = -1;
        for (size_t i = 0; i < polygon->degree(); i++)
        {
            if (polygon->edgeId(i) == seed.edgeId)
            {
                edgeIndexInFace = i;
                break;
            }
        }
        assert(edgeIndexInFace != -1);
        MyVertex::Index vIdInPlane = polygon->vertexId((edgeIndexInFace + 2) % polygon->degree());

        // find a bounding plane
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
                    meshSets.insert(neigh.neighborEdgeId);
                }
            }
            edge.neighbor->swap(newNeighbor);
        }

        XPlane boundPlane;
        bool flag = false;
        for (auto &neigh : *edge.neighbor)
        {
            if (neigh.type == NeighborInfo::Edge)
            {
                for (auto fItr = MyEdge::FaceIterator(xedge(neigh.neighborEdgeId));
                    fItr; fItr.incrementToTriangle())
                {
                    boundPlane = ((Triangle*)fItr.face())->supportingPlane();
                    if (orientation(boundPlane, xvertex(vIdInPlane)) == ON_ORIENTED_BOUNDARY)
                    {
                        flag = true;
                        break;
                    }
                }
            }
            else
            {
                assert(neigh.type == NeighborInfo::Face);
                boundPlane = neigh.pTrangle->supportingPlane();
                if (orientation(boundPlane, xvertex(vIdInPlane)) == ON_ORIENTED_BOUNDARY)
                    break;
            }
            if (flag) break;
        }

        assert(boundPlane.isValid());

        // correct the direction of bounding plane
        XLine edgeLine(polygon->supportingPlane(), boundPlane);
        if (linearOrder(edgeLine, xvertex(polygon->vertexId(edgeIndexInFace + 1) % polygon->degree()),
            xvertex(polygon->vertexId(edgeIndexInFace))) < 0)
            boundPlane.inverse();

        // pick a correct rep vertex
        MyVertex::Index repVertexId;
        for (size_t i = 0; i < polygon->degree(); i++)
        {
            if (orientation(boundPlane, xvertex(polygon->vertexId(i))) == ON_POSITIVE_SIDE)
            {
                repVertexId = polygon->vertexId(i);
                break;
            }
        }
        MyVertex& repVertex = xvertex(repVertexId);

        // copy the relation
        const size_t nMesh = ((FullIndicatorVector*)seed.eIndicators.get())->getNumber();
        for (size_t i = 0; i < nMesh; i++)
            relTab[i] = (Relation)seed.eIndicators.get()->at(i);

        // correct the relation
        for (auto &neigh: *edge.neighbor)
        {
            Oriented_side side;

            if (neigh.type == NeighborInfo::Edge)
            {
                std::vector<IPolygon*> faces;
                auto fItr = MyEdge::FaceIterator(xedge(neigh.neighborEdgeId));
                for (; fItr; fItr.incrementToTriangle())
                    faces.push_back(fItr.face());

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
        Relation relation;
        TestTree dummyForest;
        std::vector<Relation> relTab(nMesh);
        std::vector<MyEdge::Index> edges;

        interQueue.push(tmpSeed);
        tmpSeed.pFace->mark = SEEDED0;
        while (!interQueue.empty())
        {
            SSeed curSeed = interQueue.front();
            interQueue.pop();

            if (curSeed.pFace->mark == VISITED)
                continue;

            curMeshId = tmpSeed.pFace->meshId();
            assert(intraQueue.empty());
            curSeed.pFace->mark = SEEDED1;
            intraQueue.push(curSeed);

            while (!intraQueue.empty())
            {
                curSeed = intraQueue.front();
                intraQueue.pop();

                if (curSeed.pFace->mark == VISITED)
                    continue;

                CSGTreeNode* tree0 = copy2(tree->pRoot, curTreeLeaves);
                calcFaceIndicator(curSeed, relTab);

                relation = ParsingCSGTree(meshList[curMeshId], relTab.data(), 
                    nMesh, tree0, curTreeLeaves, dummyForest);
                assert(relation != REL_NOT_AVAILABLE || relation != REL_UNKNOWN);

                inverse = meshList[curMeshId]->inverse();
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
                    {
                        result->faces().push_back(curFace);
                        result->inverseMap.push_back(inverse);
                    }

                    edges.clear();
                    curFace->getEdges(edges);
                    assert(edges.size() == curFace->degree());
                    for (int i = 0; i < curFace->degree(); i++)
                    {
                        MyEdge& curEdge = xedge(edges[i]);
                        MyEdge::FaceIterator fItr(curEdge);
                        if (curEdge.neighbor)
                        {
                            assert(curEdge.faceCount() >= 4);
                            for (; fItr; ++fItr)
                            {
                                tmpSeed.edgeId = edges[i];
                                tmpSeed.pFace = fItr.face();
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
                                if (fItr.face()->mark < SEEDED2)
                                {
                                    faceQueue.push(fItr.face());
                                    fItr.face()->mark = SEEDED2;
                                }
                            }
                        }
                    }
                }
            }
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