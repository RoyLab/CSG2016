#include "precompile.h"
#include <list>
#include <map>
#include <set>
#include <xlogger.h>
#include "xmemory.h"
#include "RegularMesh.h"
#include "intersection.h"
#include "XStruct.hpp"
#include "XException.hpp"

namespace Boolean
{
    namespace
    {
        struct PlaneVertex
        {
            XPlane plane;
            MyVertex::Index id;
        };

        // 用于传到sortObj里面去的辅助结构
        struct EdgeAuxiliaryStructure
        {
            MyVertex::Index start, end;
            XLine line;
        };

        struct LinOrderObj
        {
            bool operator()(const PlaneVertex& a, const PlaneVertex& b) const
            {
                int type = 0;
                if (xvertex(a.id).isPlaneRep()) type += 1;
                if (xvertex(b.id).isPlaneRep()) type += 2;

                switch (type)
                {
                case 0:
                    return line.linearOrder(xvertex(a.id).point(),
                        xvertex(b.id).point()) > 0;
                case 1:
                    return a.plane.orientation(xvertex(b.id).point()) == ON_POSITIVE_SIDE;
                case 2:
                    return b.plane.orientation(xvertex(a.id).point()) == ON_NEGATIVE_SIDE;
                case 3:
                    return line.linearOrder(a.plane, b.plane) > 0;
                default:
                    throw std::exception();
                }
            }

            XLine line;
        };

        class TessGraph
        {
            struct Edge;
            struct Node;
            typedef uint32_t EdgeIndex;
            typedef uint32_t VertexIndex;
            enum Direction { D_NAN = 0, D_SEQ = 1, D_INV = 2, D_NODIR = 3 };

            struct Node
            {
                enum { INVALID_INDEX = INVALID_UINT32 };
                XR::RecursiveVector<EdgeIndex> edges;

                EdgeIndex findnext(EdgeIndex) const;
            };

            typedef std::map<uint32_t, Node> NodeMap;

            struct Edge
            {
                NodeMap::iterator v[2];
                Direction dir;
                XPlane prep;
                ExternPtr PBIRep* pbi = nullptr;

                VertexIndex startVertex() const { return v[0]->first; }
                VertexIndex endVertex() const { return v[1]->first; }

                bool checkPlaneOrientation(const Triangle*);
            };
        public:
            struct SortObject /// clockwise
            {
                bool operator() (EdgeIndex i, EdgeIndex j);
                const TessGraph* pTG;
                VertexIndex node;
            };

        public:
            TessGraph(const Triangle*);
            bool tessellate();

        protected:
            NodeMap m_nodes;
            std::vector<Edge> m_edges;
            const Triangle* mp_tri;
        };

        TessGraph::TessGraph(const Triangle *tri)
        {
            mp_tri = tri;
            Node node; Edge edge;
            for (int i = 0; i < 3; i++)
            {
                uint32_t vId = tri->vertexId(i);
                m_nodes[vId] = node;

                auto& e = tri->edge(i);
                if (!e.inscts) continue;
                for (auto vId : e.inscts->points)
                    m_nodes[vId] = node;
            }

            if (tri->inscts)
            {
                for (auto v : tri->inscts->points)
                    m_nodes[v.vId] = node;
            }

            edge.dir = D_SEQ;
            for (int i = 0; i < 3; i++)
            {
                auto& e = tri->edge(i);
                edge.prep = tri->boundingPlane(i).opposite();

                int ires = e.faceOrientation(tri);
                assert(ires != 0);
                int i0 = 0, i1 = 1;
                if (ires < 0) std::swap(i0, i1);

                if (!e.inscts)
                {
                    auto findres0 = m_nodes.find(e.ends[0]);
                    auto findres1 = m_nodes.find(e.ends[1]);

                    if (findres0 != m_nodes.end() && findres1 != m_nodes.end())
                    {
                        edge.v[i0] = findres0;
                        edge.v[i1] = findres1;
                        edge.pbi = nullptr;
                        m_edges.push_back(edge);
                        assert(edge.checkPlaneOrientation(tri));

                        edge.v[i0]->second.edges.push_back(m_edges.size() - 1);
                        edge.v[i1]->second.edges.push_back(m_edges.size() - 1);
                    }
                    else XLOG_ERROR << "Cannot find vertex on tessellating triangle " 
                        << tri->meshId() << '/' << tri->id();
                    continue;
                }

                for (auto &set : e.inscts->inscts)
                {
                    for (auto& ePBI : set.second)
                    {
                        auto findres0 = m_nodes.find(ePBI.ends[0]);
                        auto findres1 = m_nodes.find(ePBI.ends[1]);

                        if (findres0 != m_nodes.end() && findres1 != m_nodes.end())
                        {
                            edge.v[i0] = findres0;
                            edge.v[i1] = findres1;
                            assert(edge.checkPlaneOrientation(tri));

                            edge.pbi = &ePBI;
                            m_edges.push_back(edge);

                            edge.v[i0]->second.edges.push_back(m_edges.size() - 1);
                            edge.v[i1]->second.edges.push_back(m_edges.size() - 1);
                        }
                        else XLOG_ERROR << "*Cannot find vertex on tessellating triangle "
                            << tri->meshId() << '/' << tri->id();
                    }
                }
            }

            if (!tri->inscts) return;

            edge.dir = D_NODIR;
            for (auto &set : tri->inscts->inscts)
            {
                for (auto& fPBI : set.second)
                {
                    auto findres0 = m_nodes.find(fPBI.ends[0]);
                    auto findres1 = m_nodes.find(fPBI.ends[1]);

                    if (findres0 != m_nodes.end() && findres1 != m_nodes.end())
                    {
                        edge.v[0] = findres0;
                        edge.v[1] = findres1;

                        edge.prep = fPBI.vertPlane;
                        assert(edge.checkPlaneOrientation(tri));

                        edge.pbi = &fPBI;
                        m_edges.push_back(edge);
                        edge.v[0]->second.edges.push_back(m_edges.size() - 1);
                        edge.v[1]->second.edges.push_back(m_edges.size() - 1);
                    }
                    else XLOG_ERROR << "*Cannot find vertex on tessellating triangle "
                        << tri->meshId() << '/' << tri->id();
                }
            }
        }

        bool checkSeq(std::vector<uint32_t>& vec, TessGraph::SortObject& sortObj)
        {
            for (int i = 1; i < vec.size(); i++)
            {
                if (!sortObj(vec[i - 1], vec[i]) && sortObj(vec[i], vec[i - 1]))
                    return false;
            }
            return true;
        }

        bool TessGraph::tessellate()
        {
            bool error;
            if (m_edges.size() < 3 || m_nodes.size() < 3)
            {
                XLOG_ERROR << "Invalid tess graph, triangle" << mp_tri->meshId() << "/" << mp_tri->id();
                return false;
            }

            // sort all node by circular order
            SortObject sortObj = { this };
            for (auto& nPair : m_nodes)
            {
                auto &node = nPair.second;
                sortObj.node = nPair.first;
                std::quicksort(node.edges.begin(), node.edges.end(), sortObj);
                assert(checkSeq(node.edges, sortObj));
            }

            // tessellate
            typedef uint32_t VeretxIndex;
            NodeMap::iterator chead, curV;
            EdgeIndex cedge;
            std::stack<EdgeIndex> edgeStack;
            std::vector<uint32_t> loop;
            std::vector<EdgeIndex> loopEdge;
            edgeStack.push(0);
            auto pMem = MemoryManager::getInstance();

            while (!edgeStack.empty())
            {
                if (m_edges[edgeStack.top()].dir == D_NAN)
                {
                    edgeStack.pop();
                    continue;
                }

                cedge = edgeStack.top();
                edgeStack.pop();

                auto& eRef = m_edges[cedge];
                assert(eRef.dir != D_NODIR);

                loop.clear();
                loopEdge.clear();
                error = false;

                chead = (eRef.dir == D_SEQ) ? eRef.v[0] : eRef.v[1];
                curV = (eRef.dir == D_SEQ) ? eRef.v[1] : eRef.v[0];

                loop.push_back(chead->first);
                loopEdge.push_back(cedge);

                while (curV != chead)
                {
                    cedge = curV->second.findnext(cedge);
                    Edge& newEdge = m_edges[cedge];

                    Direction dir = D_SEQ;
                    if (newEdge.v[0] != curV)
                        dir = D_INV;

                    assert(newEdge.dir > 0 && (newEdge.dir & dir));
                    if (newEdge.dir < 1)
                    {
                        XLOG_ERROR << "Tessellation error: " << mp_tri->meshId() << "/" << mp_tri->id();
                        error = true;
                        break;
                    }
                    newEdge.dir = Direction(newEdge.dir - dir);
                    loop.push_back(curV->first);
                    loopEdge.push_back(cedge);

                    auto& eRef = m_edges[cedge];
                    curV = (dir == D_SEQ) ? newEdge.v[1] : newEdge.v[0];

                    if (newEdge.dir != D_NAN)
                        edgeStack.push(cedge);
                }

                if (error) continue;

                // add subpolygon
                const uint32_t n = loop.size();
                SubPolygon *spoly = new SubPolygon(mp_tri->meshId(), n);
                spoly->constructFromVertexList(loop.begin(), loop.end());

                // add neighborInfo
                for (uint32_t i = 0; i < n; i++)
                {
                    auto* pbi = m_edges[loopEdge[i]].pbi;
                    if (pbi && !pbi->neighbor.empty())
                    {
                        auto &nSourse = pbi->neighbor;
                        auto &nTarget = spoly->edge(i).neighbor;
                        if (!nTarget)
                            nTarget = new std::vector<NeighborInfo>;
                        nTarget->insert(nTarget->end(), nSourse.begin(), nSourse.end());
                    }
                }

                spoly->sPlane = mp_tri->supportingPlane();
                pMem->addSubPolygon(spoly);
            }
            return true;
        }

        bool TessGraph::SortObject::operator()(EdgeIndex i, EdgeIndex j)
        {
            if (i == j) return false;

            auto &ei = pTG->m_edges[i];
            auto &ej = pTG->m_edges[j];

            auto sp = pTG->mp_tri->supportingPlane();
            XPlane ep[2] = { ei.prep, ej.prep };

            if (ei.startVertex() != node)
                ep[0].inverse();

            if (ej.startVertex() != node)
                ep[1].inverse();

            Real res = sign(sp, ep[0], ep[1]);
            //assert(res != Real(0));
            return res < Real(0);
        }

        TessGraph::EdgeIndex TessGraph::Node::findnext(EdgeIndex now) const
        {
            for (int i = 0; i < edges.size(); i++)
            {
                if (edges[i] == now)
                    return edges[i + 1];
            }
            return INVALID_INDEX;
        }
        bool TessGraph::Edge::checkPlaneOrientation(const Triangle *pTri)
        {
            XLine line(pTri->supportingPlane(), prep);
            if (linearOrder(line, xvertex(v[0]->first),
                xvertex(v[1]->first)) > 0)
                return true;

            return false;
        }

        struct FacePBITessData
        {
            decltype(FaceInsctData::inscts)::mapped_type::iterator ptr;
            decltype(FaceInsctData::inscts)::mapped_type* pContainer = nullptr;

            std::vector<PlaneVertex> points;
        };
    }

    void EdgeInsctData::refine(void* pData)
    {
        if (isRefined()) return;

        std::vector<PlaneVertex> seqs(points.size());
        auto vItr = points.begin();
        for (int i = 0; i < points.size(); i++, vItr++)
            seqs[i].id = *vItr;

        std::map<MyVertex::Index, XPlane> v2p;
        for (auto& set : inscts)
        {
            for (auto &pbi : set.second)
            {
                v2p[pbi.ends[0]] = pbi.pends[0];
                v2p[pbi.ends[1]] = pbi.pends[1];
            }
        }

        EdgeAuxiliaryStructure data = *(EdgeAuxiliaryStructure*)(pData);
        // assign each vertex a plane and the corresponding id
        for (int i = 0; i < points.size(); i++)
        {
            auto res0 = v2p.find(seqs[i].id);
            if (res0 == v2p.end())
            {
                // some vertex cannot be assigned a valid plane
                // therefore we should find it manually from Plane Triples
                auto &vRef = xvertex(seqs[i].id);
                if (vRef.isPlaneRep())
                {
                    auto& xpointRef = vRef.ppoint();
                    for (int j = 0; j < 3; j++)
                    {
                        Real fres = data.line.dot(xpointRef.plane(j));
                        if (fres == Real(0)) continue;

                        if (fres > 0)
                            seqs[i].plane = xpointRef.plane(j);
                        else if (fres < 0)
                            seqs[i].plane = xpointRef.plane(j).opposite();
                        break;
                    }
                }
            }
            else seqs[i].plane = res0->second;
        }

        LinOrderObj orderObj = { data.line };
        std::sort(seqs.begin(), seqs.end(), orderObj);
        std::vector<EdgePBI> newPbi(points.size() + 1);
        std::map<MyVertex::Index, uint32_t> idmap;
        idmap[data.start] = 0;
        idmap[data.end] = points.size() + 1;
        for (int i = 0; i < seqs.size(); i++)
        {
            idmap[seqs[i].id] = i + 1;
            newPbi[i].ends[1] = seqs[i].id;
            newPbi[i + 1].ends[0] = seqs[i].id;
            newPbi[i].pends[1] = seqs[i].plane;
            newPbi[i + 1].pends[0] = seqs[i].plane;
        }
        newPbi[0].ends[0] = data.start;
        newPbi[points.size()].ends[1] = data.end;

        // add pbi intersections' neighborInfo into newly constructed pbi.
        for (auto& set : inscts)
        {
            for (auto &pbi : set.second)
            {
                assert(idmap.find(pbi.ends[0]) != idmap.end());
                assert(idmap.find(pbi.ends[1]) != idmap.end());

                uint32_t start = idmap[pbi.ends[0]];
                uint32_t last = idmap[pbi.ends[1]];

                assert(start < last);
                for (int i = start; i < last; i++)
                {
                    assert(pbi.neighbor.size() == 1);
                    newPbi[i].neighbor.push_back(*pbi.neighbor.begin());
                }
            }
        }
        inscts.clear();
        auto& slot = inscts[INVALID_UINT32];
        for (int i = 0; i < newPbi.size(); i++)
            slot.push_back(newPbi[i]);

        bRefined = true;
    }

    bool checkIsolatedPart(Triangle* pTri, std::vector<MyVertex::Index>& strayVertices)
    {
        struct ColorVertex
        {
            std::set<MyVertex::Index> neighbors;
            int color = -1;
        } defaultCV;

        FaceInsctData* thiz = pTri->inscts;
        // initialize the graph
        std::map<MyVertex::Index, ColorVertex> data;
        for (auto& v : thiz->points)
            data[v.vId] = defaultCV;

        for (auto& pbiSet : thiz->inscts)
        {
            for (FacePBI& fPbi : pbiSet.second)
            {
                data[fPbi.ends[0]].neighbors.insert(fPbi.ends[1]);
                data[fPbi.ends[1]].neighbors.insert(fPbi.ends[0]);
            }
        }

        // init vertex set
        std::set<MyVertex::Index> edgeVertices;
        for (int i = 0; i < 3; i++)
        {
            edgeVertices.insert(pTri->vertexId(i));
            MyEdge& edge = pTri->edge(i);
            if (!edge.inscts) continue;
            for (MyVertex::Index vId : edge.inscts->points)
                edgeVertices.insert(vId);
        }

        // colorization
        int color = 0;
        uint32_t isoCount = MAX_MESH_COUNT;
        bool isIsolated = false, boolRes = false;
        std::stack<MyVertex::Index> vStack;
        std::vector<MyVertex::Index> history;
        for (auto& vItem : data)
        {
            if (vItem.second.color == -1)
            {
                ++color;
                isIsolated = true;
                assert(vStack.empty());
                vStack.push(vItem.first);
                while (!vStack.empty())
                {
                    MyVertex::Index curVId = vStack.top();
                    vStack.pop();
                    auto& curV = data[curVId];
                    if (curV.color != -1) continue;

                    history.push_back(curVId);
                    if (edgeVertices.find(curVId) != edgeVertices.end())
                        isIsolated = false;

                    curV.color = color;
                    for (MyVertex::Index vId : curV.neighbors)
                    {
                        if (data[vId].color == -1)
                            vStack.push(vId);
                    }
                }

                if (isIsolated)
                {
                    MyVertex::Index chooseVertex = INVALID_UINT32;
                    MyEdge::SIndex eId = -1;
                    for (MyVertex::Index vId : history)
                    {
                        for (auto &v : thiz->points)
                        {
                            if (v.vId == vId)
                            {
                                chooseVertex = vId;
                                eId = v.eId;
                                break;
                            }
                        }

                        if (eId >= 0) break;
                    }
                    assert(eId != -1);
                    if (eId >= 0)
                    {
                        FacePBI fPbi;
                        fPbi.ends[0] = pTri->vertexId(0);
                        fPbi.ends[1] = chooseVertex;

                        MyEdge& crossEdgeRef = xedge(eId);
                        fPbi.vertPlane = XPlane(pTri->vertex(0).point(),
                            xvertex(crossEdgeRef.ends[0]).point(), xvertex(crossEdgeRef.ends[1]).point());

                        XLine line(pTri->supportingPlane(), fPbi.vertPlane);
                        fPbi.pends[0] = pickPositiveVertical(line, xvertex(fPbi.ends[0]));
                        fPbi.pends[1] = pickPositiveVertical(line, xvertex(fPbi.ends[1]));
                        assert(line.dot(fPbi.pends[0]) > 0);
                        assert(line.dot(fPbi.pends[1]) > 0);
                        if (line.linearOrderNoCheck(fPbi.pends[0], fPbi.pends[1]) < 0)
                        {
                            std::swap(fPbi.ends[0], fPbi.ends[1]);
                            std::swap(fPbi.pends[0], fPbi.pends[1]);
                        }

                        assert(thiz->inscts.find(isoCount) == thiz->inscts.end());
                        thiz->inscts[isoCount++].push_back(fPbi);
                    }
                    else
                    {
                        assert(eId == -2 && history.size() == 1);
                        strayVertices.push_back(history[0]);
                    }
                    boolRes = true;
                }
                history.clear();
            }
        }

        return boolRes;
    }

    void removeOverlapPBI(FaceInsctData * thiz)
    {
        std::map<MyVertex::Index, std::set<FacePBI*>> data;
        std::vector<decltype(thiz->inscts[0].begin())> garbage;
        for (auto &pbiSet : thiz->inscts)
        {
            for (auto pbiItr = pbiSet.second.begin(); pbiItr != pbiSet.second.end(); ++pbiItr)
            {
                FacePBI& pbi = *pbiItr;
                bool found = false;
                for (auto alreadyHere : data[pbi.ends[0]]) // search in current graph
                {
                    if (alreadyHere->ends[0] == pbi.ends[1] ||
                        alreadyHere->ends[1] == pbi.ends[1]) // if has
                    {
                        for (NeighborInfo& nInfo : pbi.neighbor)
                        {
                            bool unique = true; // search if the neighInfo is already there, == std::find
                            for (NeighborInfo& nInfo2 : alreadyHere->neighbor)
                            {
                                if (nInfo.neighborMeshId == nInfo.neighborMeshId)
                                {
                                    unique = false;
                                    break;
                                }
                            }

                            if (unique)
                            {
                                alreadyHere->neighbor.push_back(nInfo);
                            }
                        }
                        found = true;
                        break;
                    }
                }

                if (found)
                {
                    garbage.push_back(pbiItr);
                    continue;
                }

                data[pbi.ends[0]].insert(&pbi);
                data[pbi.ends[1]].insert(&pbi);
            }

            for (auto& pbiItr : garbage)
                pbiSet.second.erase(pbiItr);
            garbage.clear();
        }
    }

    IndexPair makePbiIndex(const PBIRep* rep)
    {
        IndexPair res;
        MakeIndex(rep->ends, res);
        return res;
    }

    void FaceInsctData::resolveIntersection(Triangle* pTri, std::vector<MyVertex::Index>* strayVertices)
    {
        XPlane triSp = pTri->supportingPlane();
        std::map<IndexPair, std::shared_ptr<FacePBITessData>> tessData;
        for (auto setItr = inscts.begin(); setItr != inscts.end(); ++setItr)
        {
            auto setItr2 = setItr; ++setItr2;
            for (; setItr2 != inscts.end(); ++setItr2)
            {
                for (auto pbiItr = setItr->second.begin(); pbiItr != setItr->second.end(); ++pbiItr)
                {
                    FacePBI& pbi = *pbiItr;
                    for (auto pbiItr2 = setItr2->second.begin(); pbiItr2 != setItr2->second.end(); ++pbiItr2)
                    {
                        PlaneVertex slots[2][2]; // 需要添加的顶点
                        FacePBI& pbi2 = *pbiItr2;
                        Oriented_side side[2][2];
                        side[0][0] = orientation(pbi2.vertPlane, pbi.ends[0]);
                        side[0][1] = orientation(pbi2.vertPlane, pbi.ends[1]);
                        if (side[0][0] == side[0][1])
                        {
                            if (side[0][0] == ON_ORIENTED_BOUNDARY) // 共线情况
                            {
                                XLine line(triSp, pbi.vertPlane);
                                assert(line.dot(pbi2.vertPlane) == 0.);
                                Real dotRes = line.dot(pbi2.pends[0]);
                                assert(dotRes != 0.);
                                bool inverseLine = dotRes > 0 ? false : true;

                                PlaneVertex pbi2ends[2] = { { pbi2.pends[0], pbi2.ends[0] },{ pbi2.pends[1], pbi2.ends[1] } };
                                if (inverseLine)
                                {
                                    pbi2ends[0].plane.inverse();
                                    pbi2ends[1].plane.inverse();
                                    std::swap(pbi2ends[0], pbi2ends[1]);
                                }

                                // linear order 计算overlap
                                assert(line.linearOrder(pbi2ends[0].plane, pbi2ends[1].plane) > 0);
                                assert(line.linearOrderNoCheck(pbi2ends[0].plane, pbi2ends[1].plane) > 0);
                                if (line.linearOrderNoCheck(pbi2ends[0].plane, pbi.pends[1]) <= 0 ||
                                    line.linearOrderNoCheck(pbi.pends[0], pbi2ends[1].plane) <= 0)
                                    continue;

                                int compRes[2] = {
                                    line.linearOrderNoCheck(pbi.pends[0], pbi2ends[0].plane),
                                    line.linearOrderNoCheck(pbi.pends[1], pbi2ends[1].plane),
                                };

                                assert(!(compRes[0] == 0 && compRes[1] == 0));

                                if (compRes[0] > 0)
                                    slots[0][0] = pbi2ends[0];
                                else if (compRes[0] < 0)
                                {
                                    slots[1][0].id = pbi.ends[0];
                                    slots[1][0].plane = pbi.pends[0];
                                    if (inverseLine) slots[1][0].plane.inverse();
                                }

                                if (compRes[1] > 0)
                                {
                                    slots[1][1].id = pbi.ends[1];
                                    slots[1][1].plane = pbi.pends[1];
                                    if (inverseLine) slots[1][1].plane.inverse();
                                }
                                else if (compRes[1] < 0)
                                    slots[0][1] = pbi2ends[1];
                            }
                            else continue;
                        }
                        else
                        {
                            side[1][0] = orientation(pbi.vertPlane, pbi2.ends[0]);
                            side[1][1] = orientation(pbi.vertPlane, pbi2.ends[1]);

                            if (side[1][0] == side[1][1])
                            {
                                assert(side[1][0] != ON_ORIENTED_BOUNDARY);//如果是共线情况，在前面应当已经被测试到
                                continue; // 不相交
                            }

                            if (side[0][0] * side[0][1] == 0)
                            {
                                if (side[1][0] * side[1][1] == 0) // 相交于某个顶点，不用split
                                    continue;

                                int addedTarget = side[0][0] == ON_ORIENTED_BOUNDARY ? 0 : 1;
                                slots[1][0].id = pbi.ends[addedTarget];
                                slots[1][0].plane = pbi.vertPlane;
                                XLine(triSp, pbi2.vertPlane).makePositive(slots[1][0].plane);
                            }
                            else
                            {
                                assert(side[0][0] * side[0][1] == -1);
                                if (side[1][0] * side[1][1] == 0)
                                {
                                    int addedTarget = side[1][0] == ON_ORIENTED_BOUNDARY ? 0 : 1;
                                    slots[0][0].id = pbi2.ends[addedTarget];
                                    slots[0][0].plane = pbi2.vertPlane;
                                    XLine(triSp, pbi.vertPlane).makePositive(slots[0][0].plane);
                                }
                                else
                                {
                                    assert(side[1][0] * side[1][1] == -1);
                                    // new vertex
                                    XPlane thirdPlane = pbi2.vertPlane;
                                    XLine(triSp, pbi.vertPlane).makePositive(thirdPlane); // 似乎不需要，可以尝试注释这一句
                                    XPoint newPoint(triSp, pbi.vertPlane, thirdPlane);

                                    uint32_t* newPos;
                                    newPos = point(newPoint, -1);
                                    std::vector<uint32_t*> vecs;
                                    vecs.push_back(newPos);
                                    uint32_t minVal = *newPos;

                                    // 去所有的邻居看一看
                                    FacePBI* twoPbiPtr[2] = { &pbi, &pbi2 };

                                    for (int i = 0; i < 2; i++)
                                    {
                                        int i2 = (i + 1) % 2;
                                        for (NeighborInfo& neiInfo : twoPbiPtr[i]->neighbor)
                                        {
                                            if (neiInfo.type == NeighborInfo::Edge)
                                            {
                                                MyEdge& eRef = xedge(neiInfo.neighborEdgeId);
                                                assert(eRef.inscts);
                                                newPos = eRef.inscts->point(newPoint);
                                            }
                                            else
                                            {
                                                assert(neiInfo.type == NeighborInfo::Face);
                                                assert(neiInfo.pTrangle->inscts);
                                                int eId = -1;
                                                if (twoPbiPtr[i2]->neighbor.empty())
                                                    eId = -2;

                                                newPos = neiInfo.pTrangle->inscts->point(newPoint, eId);
                                            }
                                            if (*newPos < minVal)
                                                minVal = *newPos;
                                            vecs.push_back(newPos);
                                        }
                                    }
                                    
                                    // 所有的邻居都没有，那就真没有了
                                    if (minVal == INVALID_UINT32)
                                    {
                                        minVal = MemoryManager::getInstance()->insertVertex(newPoint);
                                    }
                                    for (uint32_t *pInt : vecs)
                                        *pInt = minVal;

                                    slots[0][0].id = minVal;
                                    slots[0][0].plane = pbi2.vertPlane;
                                    XLine(triSp, pbi.vertPlane).makePositive(slots[0][0].plane);

                                    slots[1][0].id = minVal;
                                    slots[1][0].plane = pbi.vertPlane;
                                    XLine(triSp, pbi2.vertPlane).makePositive(slots[1][0].plane);
                                }
                            }

                        }

                        if (slots[0][0].plane.isValid() || slots[0][1].plane.isValid())
                        {
                            IndexPair pbiIndex = makePbiIndex(&pbi);
                            auto pbiData = tessData.find(pbiIndex);
                            if (pbiData == tessData.end())
                            {
                                FacePBITessData *tessItem = new FacePBITessData;
                                tessItem->pContainer = &setItr->second;
                                tessItem->ptr = pbiItr;
                                auto insertRes = tessData.insert(decltype(tessData)::value_type(pbiIndex,
                                    std::shared_ptr<FacePBITessData>(tessItem)));

                                assert(insertRes.second);
                                pbiData = insertRes.first;
                            }

                            if (slots[0][0].plane.isValid())
                                pbiData->second->points.push_back(slots[0][0]);

                            if (slots[0][1].plane.isValid())
                                pbiData->second->points.push_back(slots[0][1]);
                        }

                        if (slots[1][0].plane.isValid() || slots[1][1].plane.isValid())
                        {
                            IndexPair pbiIndex2 = makePbiIndex(&pbi2);
                            auto pbiData = tessData.find(pbiIndex2);
                            if (pbiData == tessData.end())
                            {
                                FacePBITessData *tessItem = new FacePBITessData;
                                tessItem->pContainer = &setItr2->second;
                                tessItem->ptr = pbiItr2;
                                auto insertRes = tessData.insert(decltype(tessData)::value_type(pbiIndex2,
                                    std::shared_ptr<FacePBITessData>(tessItem)));

                                assert(insertRes.second);
                                pbiData = insertRes.first;
                            }

                            if (slots[1][0].plane.isValid())
                                pbiData->second->points.push_back(slots[1][0]);

                            if (slots[1][1].plane.isValid())
                                pbiData->second->points.push_back(slots[1][1]);
                        }
                    }
                }
            }
        }

        if (strayVertices)
        {
            for (MyVertex::Index strayV : *strayVertices)
            {
                if (!xvertex(strayV).isPlaneRep()) continue; // 不合理，但是先这样吧
                for (auto setItr = inscts.begin(); setItr != inscts.end(); ++setItr)
                {
                    for (auto pbiItr = setItr->second.begin(); pbiItr != setItr->second.end(); ++pbiItr)
                    {
                        FacePBI& fPbi = *pbiItr;
                        if (fPbi.vertPlane.has_on(xvertex(strayV).ppoint()))
                        {
                            IndexPair pbiIndex = makePbiIndex(&fPbi);
                            auto pbiData = tessData.find(pbiIndex);
                            if (pbiData == tessData.end())
                            {
                                FacePBITessData *tessItem = new FacePBITessData;
                                tessItem->pContainer = &setItr->second;
                                tessItem->ptr = pbiItr;
                                auto insertRes = tessData.insert(decltype(tessData)::value_type(pbiIndex,
                                    std::shared_ptr<FacePBITessData>(tessItem)));

                                assert(insertRes.second);
                                pbiData = insertRes.first;
                            }

                            XLine line(pTri->supportingPlane(), fPbi.vertPlane);
                            line.pickPositiveVertical(xvertex(strayV).ppoint());
                            pbiData->second->points.push_back(
                                PlaneVertex{ line.pickPositiveVertical(xvertex(strayV).ppoint()), strayV });
                            break; // 如果需要加两个以上的fpbi的话，这个点应该已经被探测出来了
                        }
                    }
                }
            }
        }

        for (auto &pPair : tessData)
        {
            FacePBITessData* pData = pPair.second.get();
            XLine line(triSp, pData->ptr->vertPlane);
            LinOrderObj orderObj = { line };

            auto &inserted = pData->points;
            std::sort(inserted.begin(), inserted.end(), orderObj);
            inserted.erase(std::unique(inserted.begin(), inserted.end(),
                [](const PlaneVertex &a, const PlaneVertex &b)->bool {
                return a.id == b.id;
            }), inserted.end());

            std::vector<FacePBI> newPbi(inserted.size() + 1, *pData->ptr);
            std::map<MyVertex::Index, uint32_t> idmap;
            idmap[pData->ptr->ends[0]] = 0;
            idmap[pData->ptr->ends[1]] = inserted.size() + 1;
            for (int i = 0; i < inserted.size(); i++)
            {
                idmap[inserted[i].id] = i + 1;
                newPbi[i].ends[1] = inserted[i].id;
                newPbi[i + 1].ends[0] = inserted[i].id;
                newPbi[i].pends[1] = inserted[i].plane;
                newPbi[i + 1].pends[0] = inserted[i].plane;
            }
            newPbi[0].ends[0] = pData->ptr->ends[0];
            newPbi[inserted.size()].ends[1] = pData->ptr->ends[1];

            pData->pContainer->erase(pData->ptr); 
            pData->pContainer->insert(pData->pContainer->end(), newPbi.begin(), newPbi.end());
            newPbi.clear();
        }

        removeOverlapPBI(this);
    }


    void FaceInsctData::refine(void* pData)
    {
        removeOverlapPBI(this);
        Triangle* pTri = reinterpret_cast<Triangle*>(pData);
        if (!isRefined() && inscts.size() >= 2)
            resolveIntersection(pTri);

        std::vector<MyVertex::Index> strayVertices;
        if (checkIsolatedPart(pTri, strayVertices))
            resolveIntersection(pTri, &strayVertices);

        bRefined = true;
    }

    void tessellation(std::vector<RegularMesh*>& meshes)
    {
        auto pMem = MemoryManager::getInstance();
        auto& inscts = intersectTriangles();
        for (Triangle* pTri : inscts)
        {
            assert(pTri->isAdded4Tess());
            if (pTri->inscts)
                pTri->inscts->refine((void*)pTri);

            for (int i = 0; i < 3; i++)
            {
                EdgeAuxiliaryStructure data = { pTri->edge(i).ends[0] , pTri->edge(i).ends[1] };
                if (pTri->edge(i).faceOrientation(pTri) > 0)
                    data.line = XLine(pTri->supportingPlane(), pTri->boundingPlane(i).opposite());
                else
                    data.line = XLine(pTri->supportingPlane(), pTri->boundingPlane(i));

                if (pTri->edge(i).inscts)
                    pTri->edge(i).inscts->refine((void*)&data);
            }

            TessGraph tg(pTri);
            if (tg.tessellate())
                pTri->invalidate();
        }
    }
}