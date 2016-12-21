#include "precompile.h"
#include <map>
#include <set>
#include <stack>

#include <xlogger.h>
#include <XStruct.hpp>

#include "hybrid_geometry.h"
#include "xmemory.h"
#include "intersection.h"

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

namespace Boolean
{
    namespace
    {
        class TessGraph
        {
            friend class ResolveObj;
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
                int component_index; // default -1
                VertexIndex vertex_id;
                EdgeSIndex coincident_edge_id;
                XR::RecursiveVector<ConnectionIndex> connections;
            };

            struct Connection
            {
                int                 component_index; // default -1
                NodeIndex           ends[2];
                Direction           dir;
                XPlane              prep;
                ExternPtr PbiRep*   pbi; // for neighborhood assgin

                int                 loop_index_left;
                int                 loop_index_right;
            };

            struct IntersectionResult
            {
                bool is_con_type;
                union
                {
                    ConnectionIndex con_idx;
                    NodeIndex node_idx;
                };
                XPlane plane;
            };

            typedef std::vector<ConnectionIndex> ConnectionLoop;
            typedef std::vector<NodeIndex> NodeLoop;

            struct Loop
            {
                ConnectionLoop cloop;
                NodeLoop nloop;
                bool nested;
            };
            typedef std::vector<Loop> Component;

            struct LoopLocation
            {
                int component_idx;
                int loop_idx;
            };

            class OuterInnerTable
            {
            public:
                OuterInnerTable(size_t n);
                void set_sibling_or_outer(LoopLocation outer, LoopLocation inner);
                bool is_inner_or_sibling(int inner, int outer) const;
                int get_outer() const;

                class FaceIterator
                {
                public:
                    FaceIterator(OuterInnerTable*);
                    void operator++();
                    operator bool() const;
                    void retrieve(std::vector<Loop*>&);
                };

                FaceIterator face_begin() { return FaceIterator(this); }
            };

        private:
            ConnectionIndex theFirstConnection() const { return 0; }
            //VertexIndex end_point(const Connection& con, int i) const;
            ConnectionIndex findnext(const Node& node, ConnectionIndex cur) const;
            bool checkPlaneOrientation(const Connection& con, const XPlane& sp) const;
            void nestedResolveAndExtractFaces(std::vector<Component*>&) const;
            void extract_face(const Loop& loop) const;
            void extract_complex_face(const std::vector<Loop*>&) const;
            void resolveIsolatedVertices() const;
            NodeIndex get_vertex_node_id(int i) const { return i; }
            int get_loop_from_cross_point(const IntersectionResult&, LoopLocation& loop_loc, bool inverse = false) const;

            Loop& get_loop(std::vector<Component*>& comp, LoopLocation& loc) const
            {
                return comp[loc.component_idx]->at(loc.loop_idx);
            }

            VertexIndex get_vertex_idx(NodeIndex node_idx) const
            {
                return nodes_[node_idx].vertex_id;
            }

            void addConnection(std::map<VertexIndex, NodeIndex>&, VertexIndex a, VertexIndex b, PbiRep* pbi, Direction dir, XPlane prep);

            std::vector<Node>       nodes_;
            std::vector<Connection> connections_;
            const Triangle*         triangle_;
            int                     face_pbi_offset_;
        };


        class ResolveObj :
            public IResolveLineIntersection<
            PlanePoint, cyPointT,
            MyVertex, MyVertex
            >
        {
        public:
            ResolveObj() : node_ids_{ -1,-1 } {}

            void resolve(Oriented_side side[4], const XPlane& a, const XPlane&b,
                const PlanePoint& a0, const cyPointT& a1,
                const MyVertex& b0, const MyVertex& b1);

            bool is_vertex_intersection() const { return insct_node_id > -1; }
            VertexIndex get_intersect_node_id() const
            {
#ifdef XR_DEBUG
                if (!is_vertex_intersection())
                {
                    throw 1;
                }
#endif
                return static_cast<VertexIndex>(insct_node_id);
            }

            void set_node_id(TessGraph::NodeIndex a, TessGraph::NodeIndex b)
            {
                node_ids_[0] = a;
                node_ids_[1] = b;
            }

        private:
            VertexSIndex insct_node_id;
            TessGraph::NodeIndex node_ids_[2];
        };

        // impl start
        TessGraph::TessGraph(const Triangle *tri)
        {
            triangle_ = tri;
            Node node = { -1 };
            std::map<VertexIndex, NodeIndex> vertex_node_dict;
            for (int i = 0; i < 3; i++)
            {
                node.vertex_id = triangle_->vertexId(i);
                nodes_.push_back(node);
                vertex_node_dict[node.vertex_id] = static_cast<int>(nodes_.size() - 1);

                MyEdge& e = triangle_->edge(i);
                if (!e.inscts) continue;
                for (auto vertex : e.inscts->points)
                {
                    node.vertex_id = vertex.vertex_idx;
                    nodes_.push_back(node);
                    vertex_node_dict[node.vertex_id] = static_cast<int>(nodes_.size() - 1);
                }
            }

            if (triangle_->inscts)
            {
                for (auto v : triangle_->inscts->points)
                {
                    node.vertex_id = v.vId;
                    node.coincident_edge_id = v.eId;
                    nodes_.push_back(node);
                    vertex_node_dict[node.vertex_id] = static_cast<int>(nodes_.size() - 1);
                }
            }

            for (int i = 0; i < 3; i++)
            {
                auto& e = triangle_->edge(i);
                XPlane prep = triangle_->boundingPlane(i).opposite();
                int ires = e.faceOrientation(triangle_);
                assert(ires != 0);
                int i0 = 0, i1 = 1;
                if (ires < 0) std::swap(i0, i1);

                if (!e.inscts)
                {
                    addConnection(vertex_node_dict, e.ends[i0], e.ends[i1], nullptr, D_SEQ, prep);
                }
                else
                {
                    auto pbi_itr = e.inscts->pbi_begin();
                    while (pbi_itr)
                    {
                        addConnection(vertex_node_dict, pbi_itr->ends[i0],
                            pbi_itr->ends[i1], pbi_itr.pointer(), D_SEQ, prep);
                    }
                }
            }

            if (!triangle_->inscts) return;

            auto pbi_itr = triangle_->inscts->pbi_begin();
            while (pbi_itr)
            {
                addConnection(vertex_node_dict, pbi_itr->ends[0],
                    pbi_itr->ends[1], pbi_itr.pointer(), D_NODIR, pbi_itr->vertPlane);
            }

            // sort all node by circular order
            struct CircularItem : public CircularOrderObjDefaultItem
            {
                int id;
            };

            std::vector<CircularItem> sort_array;
            CircularOrderObj<CircularItem> sortObj(triangle_->supportingPlane());
            for (int node_idx = 0; node_idx < nodes_.size(); node_idx++)
            {
                Node& node = nodes_[node_idx];
                if (node.connections.size() <= 2) continue;

                // construct sort_array
                sort_array.resize(node.connections.size());
                for (int j = 0; j < node.connections.size(); j++)
                {
                    sort_array[j].id = j;
                    Connection& con_ref = connections_[node.connections[j]];
                    sort_array[j].prep = con_ref.prep;
                    if (con_ref.ends[0] != node_idx)
                    {
                        sort_array[j].prep.inverse();
                    }
                }

                // sort
                std::quicksort(sort_array.begin(), sort_array.end(), sortObj);
                assert(sortObj.checkSeq(sort_array));

                // assign back
                std::vector<ConnectionIndex> new_cons(node.connections.size());
                for (int j = 0; j < node.connections.size(); j++)
                {
                    new_cons[j] = node.connections[sort_array[j].id];
                }
                node.connections.swap(new_cons);
            }
        }

        bool TessGraph::tessellate()
        {
            std::stack<ConnectionIndex> edge_stack;
            std::vector<Component*> components;
            const int con_sz = static_cast<int>(connections_.size());
            int component_index = -1;
            for (int i = 0; i < con_sz; i++)
            {
                if (connections_[i].dir == D_NAN) continue;
                edge_stack.push(i);

                Component* cur_component = new Component();
                components.push_back(cur_component);
                ++component_index;
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
                    cur_loop.nested = false;
                    int loop_index = static_cast<int>(cur_component->size());
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
                        nodes_[cur_node_idx].component_index = component_index;
                        cur_con.component_index = component_index;
                        if (sub_dir == D_SEQ)
                        {
                            cur_con.loop_index_left = loop_index;
                        }
                        else
                        {
                            cur_con.loop_index_right = loop_index;
                        }

                        // end loop
                        if (cur_node_idx == headnode) break;

                        // update
                        cur_con_idx = findnext(nodes_[cur_node_idx], cur_con_idx);
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
                for (Loop& loop : *components[0])
                {
                    extract_face(loop);
                }
                return;
            }

            Component* outer_component = components[0];
            MyVertex& anchor_vertex = triangle_->vertex(0);

            //std::vector<std::pair<IntersectionResult, IntersectionResult>>
            //    insct_results(components.size() - 1);

            OuterInnerTable rel_table(components.size());
            std::vector<int> nested_relation(components.size(), -1);
            std::vector<std::set<int>> nested_faces;
            for (int i = 1; i < components.size(); i++)
            {
                Component* cur_component = components[i];

                // choose a valid point and generate a test pbi
                NodeIndex chosen_node_idx = -1;
                for (Loop& loop : *cur_component)
                {
                    for (NodeIndex node_idx : loop.nloop)
                    {
                        if (nodes_[node_idx].coincident_edge_id >= 0)
                        {
                            chosen_node_idx = node_idx;
                            break;
                        };
                    }
                    if (chosen_node_idx >= 0) break;
                }

                assert(chosen_node_idx != -1);
                EdgeIndex chosen_edge_idx = nodes_[chosen_node_idx].coincident_edge_id;
                MyEdge& chosen_edge = xedge(chosen_edge_idx);
                VertexIndex chosen_vertex_idx = nodes_[chosen_node_idx].vertex_id;
                MyVertex& chosen_vertex = xvertex(chosen_vertex_idx);

                XPlane split_plane(
                    XPlaneBase(anchor_vertex.point(),
                        xvertex(chosen_edge.ends[0]).point(),
                        xvertex(chosen_edge.ends[1]).point()
                    )
                );

                // reverse split_plane if necessary
                PlaneLine split_line(triangle_->supportingPlane(), split_plane);
                if (split_line.linear_order(chosen_vertex.ppoint(), anchor_vertex.point()) < 0)
                {
                    split_plane.inverse();
                    split_line.inverse();
                }

                // intersection test init
                std::pair<IntersectionResult, IntersectionResult> insct_pair;

                insct_pair.first.plane = split_line.
                    pick_positive_vertical_plane(chosen_vertex.ppoint());

                insct_pair.first.is_con_type = false;
                insct_pair.first.node_idx = chosen_node_idx;

                insct_pair.second.plane = split_line.
                    pick_positive_vertical_plane(anchor_vertex.point());

                insct_pair.second.is_con_type = false;
                insct_pair.second.node_idx = get_vertex_node_id(0);

                int cur_component_idx = nodes_[chosen_node_idx].component_index;

                // decide the farest cross point with the same color
                // to determine the outerloop of the current component
                for (int i = face_pbi_offset_; i < connections_.size(); i++)
                {
                    assert(connections_[i].pbi &&
                        connections_[i].pbi->is_derived_cls());

                    const Connection& test_con = connections_[i];
                    if (test_con.component_index != cur_component_idx) continue;

                    FacePbi* pbi_itr = reinterpret_cast<FacePbi*>(test_con.pbi);

                    ResolveObj resolve_obj;
                    resolve_obj.set_node_id(test_con.ends[0], test_con.ends[1]);
                    LineInsctResultType insct_res = plane_based_line_intersection(
                        split_plane, pbi_itr->vertPlane,
                        chosen_vertex.ppoint(), anchor_vertex.point(),
                        xvertex(pbi_itr->ends[0]), xvertex(pbi_itr->ends[1]),
                        &resolve_obj
                    );

                    if (insct_res == Intersect)
                    {
                        if (!resolve_obj.is_vertex_intersection()) // pbi insct
                        {
                            XPlane test_plane = pbi_itr->vertPlane;
                            split_line.make_positive(test_plane);
                            if (split_line.linear_order_unsafe(
                                insct_pair.first.plane, test_plane) > 0)
                            {
                                insct_pair.first.is_con_type = true;
                                insct_pair.first.con_idx = i;
                                insct_pair.first.plane = test_plane;
                            }
                        }
                        else // vertex insct
                        {
                            MyVertex& test_point = xvertex(nodes_[resolve_obj.get_intersect_node_id()].vertex_id);
                            XPlane test_plane = split_line.pick_positive_vertical_plane(test_point.ppoint());
                            if (split_line.linear_order_unsafe(
                                insct_pair.first.plane, test_plane) > 0)
                            {
                                insct_pair.first.is_con_type = false;
                                insct_pair.first.node_idx = resolve_obj.get_intersect_node_id();
                                insct_pair.first.plane = test_plane;
                            }
                        }
                    }
                    ++pbi_itr;
                }

                // decide the closest cross point related the outerloop that 
                // does not belong to the inner or siblings as the father (however still can be siblings)
                for (int i = face_pbi_offset_; i < connections_.size(); i++)
                {
                    assert(connections_[i].pbi &&
                        connections_[i].pbi->is_derived_cls());

                    const Connection& test_con = connections_[i];
                    if (rel_table.is_inner_or_sibling(test_con.component_index, i)) continue;

                    FacePbi* pbi_itr = reinterpret_cast<FacePbi*>(test_con.pbi);

                    ResolveObj resolve_obj;
                    resolve_obj.set_node_id(test_con.ends[0], test_con.ends[1]);

                    LineInsctResultType insct_res = plane_based_line_intersection(
                        split_plane, pbi_itr->vertPlane,
                        chosen_vertex.ppoint(), anchor_vertex.point(),
                        xvertex(pbi_itr->ends[0]), xvertex(pbi_itr->ends[1]),
                        &resolve_obj
                    );

                    if (insct_res == Intersect)
                    {
                        if (!resolve_obj.get_intersect_node_id()) // pbi insct
                        {
                            XPlane test_plane = pbi_itr->vertPlane;
                            split_line.make_positive(test_plane);
                            if (split_line.linear_order_unsafe(
                                insct_pair.first.plane, test_plane) > 0 &&
                                split_line.linear_order_unsafe(
                                    insct_pair.second.plane, test_plane) < 0
                                )
                            {
                                insct_pair.second.is_con_type = true;
                                insct_pair.second.con_idx = i;
                                insct_pair.second.plane = test_plane;
                            }
                        }
                        else // vertex insct
                        {
                            MyVertex& test_point = xvertex(nodes_[resolve_obj.get_intersect_node_id()].vertex_id);
                            XPlane test_plane = split_line.pick_positive_vertical_plane(test_point.ppoint());
                            if (split_line.linear_order_unsafe(
                                insct_pair.first.plane, test_plane) > 0 &&
                                split_line.linear_order_unsafe(
                                    insct_pair.second.plane, test_plane) < 0
                                )
                            {
                                insct_pair.second.is_con_type = false;
                                insct_pair.second.node_idx = resolve_obj.get_intersect_node_id();
                                insct_pair.second.plane = test_plane;
                            }
                        }
                    }
                    ++pbi_itr;
                } // loop pbi

                // we have get the connection pair, retreive the related loop
                std::pair<LoopLocation, LoopLocation> loop_loc_pair;
                get_loop_from_cross_point(insct_pair.first, loop_loc_pair.first);
                get_loop_from_cross_point(insct_pair.second, loop_loc_pair.second, true);

                rel_table.set_sibling_or_outer(loop_loc_pair.second, loop_loc_pair.first);

                assert(i == loop_loc_pair.first.component_idx &&
                    i != loop_loc_pair.second.component_idx);

                Loop& loop_first = get_loop(components, loop_loc_pair.first);
                Loop& loop_second = get_loop(components, loop_loc_pair.second);
                loop_first.nested = loop_second.nested = true;

            } // loop components

            for (int i = 0; i < components.size(); i++)
            {
                for (Loop& loop : *components[i])
                {
                    if (!loop.nested)
                    {
                        extract_face(loop);
                    }
                }
            }

            std::vector<Loop*> complex_faces;
            auto face_itr = rel_table.face_begin();
            while (face_itr)
            {
                face_itr.retrieve(complex_faces);
                extract_complex_face(complex_faces);
                complex_faces.clear();
            }
        }

        void TessGraph::extract_face(const Loop& loop) const
        {
            // add subpolygon
            const int degree = static_cast<int>(loop.cloop.size());
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
                PbiRep* pbi = connections_[loop.cloop[i]].pbi;
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

        void TessGraph::resolveIsolatedVertices() const
        {
            bool first_detect = true;
            for (int i = 0; i < nodes_.size(); i++)
            {
                if (!(nodes_[i].component_index == -1))
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

        void TessGraph::addConnection(std::map<VertexIndex, NodeIndex>& dict, VertexIndex a, VertexIndex b, PbiRep * pbi, Direction dir, XPlane prep)
        {
            auto findres0 = dict.find(a);
            auto findres1 = dict.find(b);
            assert(findres0 != dict.end() && findres1 != dict.end());

            Connection con = { -1 };
            con.ends[0] = findres0->second;
            con.ends[1] = findres1->second;
            con.pbi = pbi;
            con.prep = prep;
            connections_.push_back(std::move(con));
            assert(checkPlaneOrientation(con, triangle_->supportingPlane()));

            ConnectionIndex con_idx = static_cast<int>(connections_.size() - 1);
            nodes_[con.ends[0]].connections.push_back(con_idx);
            nodes_[con.ends[1]].connections.push_back(con_idx);
        }


        TessGraph::ConnectionIndex TessGraph::findnext(const Node& node, ConnectionIndex cur) const
        {
            for (int i = 0; i < node.connections.size(); i++)
            {
                if (node.connections[i] == cur)
                    return node.connections[i + 1];
            }
            return -1;
        }

        bool TessGraph::checkPlaneOrientation(const Connection& con, const XPlane& sp) const
        {
            PlaneLine line(sp, con.prep);

            if (linear_order(line,
                xvertex(get_vertex_idx(con.ends[0])),
                xvertex(get_vertex_idx(con.ends[1]))
            ) > 0)
                return true;

            return false;
        }

        void ResolveObj::resolve(Oriented_side side[4],
            const XPlane & plane_a, const XPlane & plane_b,
            const PlanePoint & a0, const cyPointT & a1,
            const MyVertex & b0, const MyVertex & b1)
        {
            if (node_ids_[0] == -1)
            {
                throw 1;
            }

            if (side[2] * side[3] == 0)
            {
                if (side[2] == ON_ORIENTED_BOUNDARY)
                {
                    // intersect on chosen vertex
                    insct_node_id = node_ids_[0];
                }
                else
                {
                    assert(side[1] == ON_ORIENTED_BOUNDARY);
                    // intersect on anchor vertex
                    insct_node_id = node_ids_[1];
                }
                return;
            }

            insct_node_id = -1;
            return;
        }

    }// namespace anonymous

    void tessellation(std::vector<RegularMesh*>& meshes, std::vector<Triangle*> insct_triangles)
    {
        for (Triangle* triangle : insct_triangles)
        {
            triangle->refine();

            TessGraph tg(triangle);
            if (tg.tessellate())
                triangle->invalidate();
        }
    }
}