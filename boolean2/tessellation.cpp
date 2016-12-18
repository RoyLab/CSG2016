#include "precompile.h"
#include <list>
#include <map>
#include <set>
#include <stack>

#include <xlogger.h>
#include <XStruct.hpp>
#include <XException.hpp>

#include "hybrid_geometry.h"
#include "xmemory.h"
#include "intersection.h"

namespace Boolean
{
    namespace
    {
        // 用于传到sortObj里面去的辅助结构
        //struct EdgeAuxiliaryStructure
        //{
        //    VertexIndex start, end;
        //    PlaneLine line;
        //};



        class TessGraph
        {
        public:
            TessGraph(const Triangle*);
            bool tessellate();

        private:
            enum Direction { 
                D_NAN = 0, D_SEQ = 1, 
                D_INV = 2, D_NODIR = 3 
            };

            typedef int ConnectionIndex;
            typedef int NodeIndex;

            struct Node
            {
                VertexIndex vertex_id;
                EdgeIndex coincident_edge_id;
                XR::RecursiveVector<ConnectionIndex> connections;
                int color; // default -1
            };

            struct Connection
            {
                NodeIndex           ends[2];
                Direction           dir;
                XPlane              prep;
                ExternPtr PBIRep*   pbi; // for neighborhood assgin
                int                 color; // default -1
            };

            typedef std::vector<ConnectionIndex> ConnectionLoop;
            typedef std::vector<NodeIndex> NodeLoop;

            struct Loop
            {
                ConnectionLoop cloop;
                NodeLoop nloop;
            };
            typedef std::vector<Loop> Component;

        private:
            ConnectionIndex theFirstConnection() const { return 0; }
            VertexIndex end_point(const Connection& con, int i) const;
            ConnectionIndex findnext(NodeIndex node, ConnectionIndex cur) const;
            bool checkPlaneOrientation(const Connection& con, const XPlane& sp) const;
            void nestedResolveAndExtractFaces(std::vector<Component*>&) const;
            void extractFaces(Component* component) const;
            void resolveIsolatedVertices() const;

            std::vector<Node>       nodes_;
            std::vector<Connection> connections_;
            const Triangle*         triangle_;
        };

        TessGraph::TessGraph(const Triangle *tri)
        {
            triangle_ = tri;
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
                        connections_.push_back(edge);
                        assert(edge.checkPlaneOrientation(tri));

                        edge.v[i0]->second.connections.push_back(connections_.size() - 1);
                        edge.v[i1]->second.connections.push_back(connections_.size() - 1);
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
                            connections_.push_back(edge);

                            edge.v[i0]->second.connections.push_back(connections_.size() - 1);
                            edge.v[i1]->second.connections.push_back(connections_.size() - 1);
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
                        connections_.push_back(edge);
                        edge.v[0]->second.connections.push_back(connections_.size() - 1);
                        edge.v[1]->second.connections.push_back(connections_.size() - 1);
                    }
                    else XLOG_ERROR << "*Cannot find vertex on tessellating triangle "
                        << tri->meshId() << '/' << tri->id();
                }
            }


            // sort all node by circular order
            CircularOrderObj sortObj = { this };
            for (auto& nPair : m_nodes)
            {
                auto &node = nPair.second;
                sortObj.node = nPair.first;
                std::quicksort(node.connections.begin(), node.connections.end(), sortObj);
                assert(checkSeq(node.connections, sortObj));
            }
        }

        bool TessGraph::tessellate()
        {
            std::stack<ConnectionIndex> edge_stack;

            std::vector<Component*> components;
            const int con_sz = connections_.size();
            int color = -1;
            for (int i = 0; i < con_sz; i++)
            {
                if (connections_[i].dir == D_NAN) continue;
                edge_stack.push(i);

                Component* cur_component = new Component();
                components.push_back(cur_component);
                ++color;
                while (!edge_stack.empty())
                {
                    ConnectionIndex cur_con_idx = edge_stack.top();
                    Connection& first_con = connections_[cur_con_idx];
                    edge_stack.pop();

                    if (first_con.dir == D_NAN) continue;

                    NodeIndex headnode = -1, cur_node_idx = -1;
                    if (first_con.dir != D_INV)
                    {
                        headnode = first_con.ends[0];
                        cur_node_idx = first_con.ends[1];
                    }
                    else
                    {
                        headnode = first_con.ends[1];
                        cur_node_idx = first_con.ends[0];
                    }

                    Loop cur_loop;
                    while (1)
                    {
                        Direction sub_dir = D_SEQ;
                        Connection& cur_con = connections_[cur_con_idx];
                        if (cur_con.ends[0] != cur_node_idx)
                            sub_dir = D_INV;

                        if (cur_con.dir < 1 || cur_con.dir & sub_dir)
                        {
                            XLOG_ERROR << "Tessellation error: " << triangle_->meshId() << "/" << triangle_->id();
                            throw 1;
                        }
                        cur_con.dir = Direction(cur_con.dir - sub_dir);

                        // add to stack
                        if (cur_con.dir != D_NAN)
                            edge_stack.push(cur_con_idx);

                        // add to loop
                        cur_loop.nloop.push_back(cur_node_idx);
                        cur_loop.cloop.push_back(cur_con_idx);

                        //colorize
                        nodes_[cur_node_idx].color = color;
                        cur_con.color = color;

                        // end loop
                        if (cur_node_idx == headnode) break;

                        // update
                        cur_con_idx = findnext(cur_node_idx, cur_con_idx);
                        cur_node_idx = (sub_dir == D_SEQ) ? cur_con.ends[1] : cur_con.ends[0];
                    }
                    cur_component->push_back(std::move(cur_loop));
                }
            }

            nestedResolveAndExtractFaces(components);
            resolveIsolatedVertices();

            for (Component* component : components)
            {
                delete component;
            }

            return true;
        }

        void TessGraph::nestedResolveAndExtractFaces(std::vector<Component*>& components) const
        {
            if (components.size() == 1)
            {
                extractFaces(components[0]);
                return;
            }

            Component* outer_component = components[0];
            MyVertex& anchor_vertex = triangle_->vertex(0);
            for (int i = 1; i < components.size(); i++)
            {
                Component* cur_component = components[i];
                NodeIndex chosen_node = -1;
                for (Loop& loop : *cur_component)
                {
                    for (NodeIndex node_idx : loop.nloop)
                    {
                        if (nodes_[node_idx].coincident_edge_id >= 0)
                        {
                            chosen_node = node_idx;
                            break;
                        };
                    }
                    if (chosen_node >= 0) break;
                }

                assert(chosen_node != -1);
                EdgeIndex chosen_edge_idx = nodes_[chosen_node].coincident_edge_id;
                MyEdge& chosen_edge = xedge(chosen_edge_idx);

                XPlane split_plane(
                    XPlaneBase( anchor_vertex.point(),
                    xvertex(chosen_edge.ends[0]).point(),
                    xvertex(chosen_edge.ends[1]).point()
                    )
                );

                // TODO yesterday
                struct ColorVertex
                {
                    std::set<VertexIndex> neighbors;
                    int color = -1;
                } defaultCV;

                FaceInsctData* thiz = pTri->inscts;
                // initialize the graph
                std::map<VertexIndex, ColorVertex> data;
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
                std::set<VertexIndex> edgeVertices;
                for (int i = 0; i < 3; i++)
                {
                    edgeVertices.insert(pTri->vertexId(i));
                    MyEdge& edge = pTri->edge(i);
                    if (!edge.inscts) continue;
                    for (VertexIndex vId : edge.inscts->points)
                        edgeVertices.insert(vId);
                }

                // colorization
                int color = 0;
                uint32_t isoCount = MAX_MESH_COUNT;
                bool isIsolated = false, boolRes = false;
                std::stack<VertexIndex> vStack;
                std::vector<VertexIndex> history;
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
                            VertexIndex curVId = vStack.top();
                            vStack.pop();
                            auto& curV = data[curVId];
                            if (curV.color != -1) continue;

                            history.push_back(curVId);
                            if (edgeVertices.find(curVId) != edgeVertices.end())
                                isIsolated = false;

                            curV.color = color;
                            for (VertexIndex vId : curV.neighbors)
                            {
                                if (data[vId].color == -1)
                                    vStack.push(vId);
                            }
                        }

                        if (isIsolated)
                        {
                            VertexIndex chooseVertex = INVALID_UINT32;
                            MyEdge::SIndex eId = -1;
                            for (VertexIndex vId : history)
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

                                PlaneLine line(pTri->supportingPlane(), fPbi.vertPlane);
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
        }

        void TessGraph::extractFaces(Component* component) const
        {
            for (Loop& loop : *component)
            {
                // add subpolygon
                const int degree = loop.cloop.size();
                SubPolygon *spoly = new SubPolygon(triangle_->meshId(), degree);

                std::vector<VertexIndex> vertices(degree);
                for (int i = 0; i < degree; i++)
                {
                    vertices[i] = nodes_[loop.nloop[i]].vertex_id;
                }
                spoly->constructFromVertexList(vertices.begin(), vertices.end());

                // add neighborInfo
                for (int i = 0; i < degree; i++)
                {
                    PBIRep* pbi = connections_[loop.cloop[i]].pbi;
                    if (pbi && !pbi->neighbor.empty())
                    {
                        auto &nSourse = pbi->neighbor;
                        auto &nTarget = spoly->edge(i).neighbor;
                        if (!nTarget)
                            nTarget = new std::vector<NeighborInfo>;
                        nTarget->insert(nTarget->end(), nSourse.begin(), nSourse.end());
                    }
                }

                spoly->sPlane = triangle_->supportingPlane();
                xmeshlist()[triangle_->meshId()]->faces().push_back(spoly);
            }
            
        }

        void TessGraph::resolveIsolatedVertices() const
        {
            bool first_detect = true;
            for (int i = 0; i < nodes_.size(); i++)
            {
                if (!nodes_[i].visited)
                {
                    if (first_detect)
                    {
                        first_detect = false;
                        XLOG_DEBUG << "Detect isolated vertex" << triangle_->meshId() << "/" << triangle_->id();
                        XLOG_DEBUG << "However, I did nothing.";
                    }
                }
            }
        }

        bool TessGraph::tessellate2()
        {
            bool error;
            if (connections_.size() < 3 || m_nodes.size() < 3)
            {
                XLOG_ERROR << "Invalid tess graph, triangle" << triangle_->meshId() << "/" << triangle_->id();
                return false;
            }

            // tessellate
            typedef uint32_t VeretxIndex;
            NodeMap::iterator head_node, curV;
            EdgeIndex cur_con_idx;
            std::stack<ConnectionIndex> edgeStack;
            std::vector<uint32_t> loop;
            std::vector<ConnectionIndex> cur_loop.cloop;
            edgeStack.push(theFirstConnection());
            GlobalData* pMem = GlobalData::getObject();

            while (!edgeStack.empty())
            {
                if (connections_[edgeStack.top()].dir == D_NAN)
                {
                    edgeStack.pop();
                    continue;
                }

                cur_con_idx = edgeStack.top();
                edgeStack.pop();

                auto& cur_con = connections_[cur_con_idx];
                assert(cur_con.dir != D_NODIR);

                loop.clear();
                cur_loop.cloop.clear();
                error = false;

                head_node = (cur_con.dir == D_SEQ) ? cur_con.v[0] : cur_con.v[1];
                curV = (cur_con.dir == D_SEQ) ? cur_con.v[1] : cur_con.v[0];

                loop.push_back(head_node->first);
                cur_loop.cloop.push_back(cur_con_idx);

                while (curV != head_node)
                {
                    cur_con_idx = curV->second.findnext(cur_con_idx);
                    Edge& cur_con = connections_[cur_con_idx];

                    Direction dir = D_SEQ;
                    if (cur_con.v[0] != curV)
                        dir = D_INV;

                    assert(cur_con.dir > 0 && (cur_con.dir & dir));
                    if (cur_con.dir < 1)
                    {
                        XLOG_ERROR << "Tessellation error: " << triangle_->meshId() << "/" << triangle_->id();
                        error = true;
                        break;
                    }
                    cur_con.dir = Direction(cur_con.dir - dir);
                    loop.push_back(curV->first);
                    cur_loop.cloop.push_back(cur_con_idx);

                    auto& cur_con = connections_[cur_con_idx];
                    curV = (dir == D_SEQ) ? cur_con.v[1] : cur_con.v[0];

                    if (cur_con.dir != D_NAN)
                        edgeStack.push(cur_con_idx);
                }

                if (error) continue;

                // add subpolygon
                const uint32_t n = loop.size();
                SubPolygon *spoly = new SubPolygon(triangle_->meshId(), n);
                spoly->constructFromVertexList(loop.begin(), loop.end());

                // add neighborInfo
                for (uint32_t i = 0; i < n; i++)
                {
                    auto* pbi = connections_[cur_loop.cloop[i]].pbi;
                    if (pbi && !pbi->neighbor.empty())
                    {
                        auto &nSourse = pbi->neighbor;
                        auto &nTarget = spoly->edge(i).neighbor;
                        if (!nTarget)
                            nTarget = new std::vector<NeighborInfo>;
                        nTarget->insert(nTarget->end(), nSourse.begin(), nSourse.end());
                    }
                }

                spoly->sPlane = triangle_->supportingPlane();
                pMem->addSubPolygon(spoly);
            }
            return true;
        }

        //bool TessGraph::SortObject::operator()(EdgeIndex i, EdgeIndex j)
        //{
        //    if (i == j) return false;

        //    auto &ei = pTG->m_edges[i];
        //    auto &ej = pTG->m_edges[j];

        //    auto sp = pTG->mp_tri->supportingPlane();
        //    XPlane ep[2] = { ei.prep, ej.prep };

        //    if (ei.startVertex() != node)
        //        ep[0].inverse();

        //    if (ej.startVertex() != node)
        //        ep[1].inverse();

        //    Real res = sign(sp, ep[0], ep[1]);
        //    //assert(res != Real(0));
        //    return res < Real(0);
        //}

        TessGraph::ConnectionIndex TessGraph::Node::findnext(ConnectionIndex now) const
        {
            for (int i = 0; i < connections.size(); i++)
            {
                if (connections[i] == now)
                    return connections[i + 1];
            }
            return INVALID_UINT32;
        }
        bool TessGraph::Connection::checkPlaneOrientation(const Triangle *pTri)
        {
            PlaneLine line(pTri->supportingPlane(), prep);
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

        std::map<VertexIndex, XPlane> v2p;
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
        std::map<VertexIndex, uint32_t> idmap;
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

    // strayVertices是孤立点，从其他的三面交点处传过来的？
    bool checkIsolatedPart(Triangle* pTri, std::vector<VertexIndex>& strayVertices)
    {
        struct ColorVertex
        {
            std::set<VertexIndex> neighbors;
            int color = -1;
        } defaultCV;

        FaceInsctData* thiz = pTri->inscts;
        // initialize the graph
        std::map<VertexIndex, ColorVertex> data;
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
        std::set<VertexIndex> edgeVertices;
        for (int i = 0; i < 3; i++)
        {
            edgeVertices.insert(pTri->vertexId(i));
            MyEdge& edge = pTri->edge(i);
            if (!edge.inscts) continue;
            for (VertexIndex vId : edge.inscts->points)
                edgeVertices.insert(vId);
        }

        // colorization
        int color = 0;
        uint32_t isoCount = MAX_MESH_COUNT;
        bool isIsolated = false, boolRes = false;
        std::stack<VertexIndex> vStack;
        std::vector<VertexIndex> history;
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
                    VertexIndex curVId = vStack.top();
                    vStack.pop();
                    auto& curV = data[curVId];
                    if (curV.color != -1) continue;

                    history.push_back(curVId);
                    if (edgeVertices.find(curVId) != edgeVertices.end())
                        isIsolated = false;

                    curV.color = color;
                    for (VertexIndex vId : curV.neighbors)
                    {
                        if (data[vId].color == -1)
                            vStack.push(vId);
                    }
                }

                if (isIsolated)
                {
                    VertexIndex chooseVertex = INVALID_UINT32;
                    MyEdge::SIndex eId = -1;
                    for (VertexIndex vId : history)
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

                        PlaneLine line(pTri->supportingPlane(), fPbi.vertPlane);
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
        std::map<VertexIndex, std::set<FacePBI*>> data;
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

    void FaceInsctData::resolveIntersection(Triangle* pTri, std::vector<VertexIndex>* strayVertices)
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
                                PlaneLine line(triSp, pbi.vertPlane);
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
                                PlaneLine(triSp, pbi2.vertPlane).makePositive(slots[1][0].plane);
                            }
                            else
                            {
                                assert(side[0][0] * side[0][1] == -1);
                                if (side[1][0] * side[1][1] == 0)
                                {
                                    int addedTarget = side[1][0] == ON_ORIENTED_BOUNDARY ? 0 : 1;
                                    slots[0][0].id = pbi2.ends[addedTarget];
                                    slots[0][0].plane = pbi2.vertPlane;
                                    PlaneLine(triSp, pbi.vertPlane).makePositive(slots[0][0].plane);
                                }
                                else
                                {
                                    assert(side[1][0] * side[1][1] == -1);
                                    // new vertex
                                    XPlane thirdPlane = pbi2.vertPlane;
                                    PlaneLine(triSp, pbi.vertPlane).makePositive(thirdPlane); // 似乎不需要，可以尝试注释这一句
                                    PlanePoint newPoint(triSp, pbi.vertPlane, thirdPlane);

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
                                                MyEdge& cur_con = xedge(neiInfo.neighborEdgeId);
                                                assert(cur_con.inscts);
                                                newPos = cur_con.inscts->point(newPoint);
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
                                        minVal = GlobalData::getObject()->insertVertex(newPoint);
                                    }
                                    for (uint32_t *pInt : vecs)
                                        *pInt = minVal;

                                    slots[0][0].id = minVal;
                                    slots[0][0].plane = pbi2.vertPlane;
                                    PlaneLine(triSp, pbi.vertPlane).makePositive(slots[0][0].plane);

                                    slots[1][0].id = minVal;
                                    slots[1][0].plane = pbi.vertPlane;
                                    PlaneLine(triSp, pbi2.vertPlane).makePositive(slots[1][0].plane);
                                }
                            }

                        }

                        if (slots[0][0].plane.is_valid() || slots[0][1].plane.is_valid())
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

                            if (slots[0][0].plane.is_valid())
                                pbiData->second->points.push_back(slots[0][0]);

                            if (slots[0][1].plane.is_valid())
                                pbiData->second->points.push_back(slots[0][1]);
                        }

                        if (slots[1][0].plane.is_valid() || slots[1][1].plane.is_valid())
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

                            if (slots[1][0].plane.is_valid())
                                pbiData->second->points.push_back(slots[1][0]);

                            if (slots[1][1].plane.is_valid())
                                pbiData->second->points.push_back(slots[1][1]);
                        }
                    }
                }
            }
        }

        if (strayVertices)
        {
            for (VertexIndex strayV : *strayVertices)
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

                            PlaneLine line(pTri->supportingPlane(), fPbi.vertPlane);
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
            PlaneLine line(triSp, pData->ptr->vertPlane);
            LinOrderObj orderObj = { line };

            auto &inserted = pData->points;
            std::sort(inserted.begin(), inserted.end(), orderObj);
            inserted.erase(std::unique(inserted.begin(), inserted.end(),
                [](const PlaneVertex &a, const PlaneVertex &b)->bool {
                return a.id == b.id;
            }), inserted.end());

            std::vector<FacePBI> newPbi(inserted.size() + 1, *pData->ptr);
            std::map<VertexIndex, uint32_t> idmap;
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

        std::vector<VertexIndex> strayVertices;
        if (checkIsolatedPart(pTri, strayVertices))
            resolveIntersection(pTri, &strayVertices);

        bRefined = true;
    }

    void tessellation(std::vector<RegularMesh*>& meshes, std::vector<Triangle*> insct_triangles)
    {
        auto pMem = GlobalData::getObject();
        auto& inscts = insct_triangles;
        for (Triangle* pTri : inscts)
        {
            assert(pTri->add_as_insct_triangle);
            if (pTri->inscts)
                pTri->inscts->refine((void*)pTri);

            for (int i = 0; i < 3; i++)
            {
                EdgeAuxiliaryStructure data = { pTri->edge(i).ends[0] , pTri->edge(i).ends[1] };
                if (pTri->edge(i).faceOrientation(pTri) > 0)
                    data.line = PlaneLine(pTri->supportingPlane(), pTri->boundingPlane(i).opposite());
                else
                    data.line = PlaneLine(pTri->supportingPlane(), pTri->boundingPlane(i));

                if (pTri->edge(i).inscts)
                    pTri->edge(i).inscts->refine((void*)&data);
            }

            TessGraph tg(pTri);
            if (tg.tessellate())
                pTri->invalidate();
        }
    }
}