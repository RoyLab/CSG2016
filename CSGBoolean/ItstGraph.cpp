#include <CGAL\Cartesian_matrix.h>
#include <queue>
#include <list>
#include "ItstGraph.h"


namespace CSG
{
    class PBPointCompare
    {
        typedef Plane_ext<K> Plane;
        typedef PBPoint<K> Point;
        COMMON_PROPERTY(Plane, p);
        COMMON_PROPERTY(Plane, q);
    public:
        PBPointCompare(Plane& sp, Plane& bp)
        {
            m_p = sp; m_q = bp;
        }

        bool operator()(Point& a, Point& b)
        {
            Plane pa = a.getPlanes()[2];
            Plane pb = b.getPlanes()[2]; //之所以是2,因为2一定是个非支撑平面
            
            CGAL::Sign s = CGAL::sign_of_determinant(m_q.a(), m_q.b(), m_q.c(),
                m_p.a(), m_p.b(), m_p.c(), pa.a(), pa.b(), pa.c());
            if (s == CGAL::NEGATIVE) pa = pa.opposite();

            s = CGAL::sign_of_determinant(m_q.a(), m_q.b(), m_q.c(),
                m_p.a(), m_p.b(), m_p.c(), pb.a(), pb.b(), pb.c());
            if (s == CGAL::NEGATIVE) pb = pb.opposite();

            double sign = orientation(m_q, m_p, pa, pb);
            if (sign > 0.0) return true;
            else return false;
        }

        bool operator()(VProxyItr& a, VProxyItr& b)
        {
            // a before than b
            auto aa = a.pointer()->pos;
            auto bb = b.pointer()->pos;

#ifdef _DEBUG
            aa.computeCoord();
            bb.computeCoord();
            CGAL::Vector_3<K> p0 = m_q.orthogonal_vector();
            CGAL::Vector_3<K> p1 = m_p.orthogonal_vector();
            auto &positive = CGAL::cross_product(p0, p1);

            CGAL::Point_3<K> posa(aa.getCoord());
            CGAL::Point_3<K> posb(bb.getCoord());

            double res2 = (posa - CGAL::ORIGIN) * positive - (posb - CGAL::ORIGIN) * positive;
            bool res = operator()(aa, bb);

            assert(res && res2 < 1e-5 || !res && res2 > -1e-5);

            return res;
#endif
        }
    };

    void ItstGraph::addEdge(Edge& e)
    {
        m_edges.push_back(e);
        size_t eId = m_edges.size() - 1;
        auto res0 = m_maps.find(e.vId[0]);
        auto res1 = m_maps.find(e.vId[1]);

        assert(res0 != m_maps.end());
        assert(res1 != m_maps.end());

        auto &node0 = m_nodes[res0->second];
        auto &node1 = m_nodes[res1->second];

        Adj adjacent;
        adjacent.eId = eId;

        adjacent.vId = m_maps[e.vId[1]];
        adjacent.startId = m_maps[e.vId[0]];
        adjacent.visited = false;
        if (!(e.direction & SEQ))
            adjacent.visited = true;
        node0.edges.push_back(adjacent);

        adjacent.vId = m_maps[e.vId[0]];
        adjacent.startId = m_maps[e.vId[1]];
        adjacent.visited = false;
        if (!(e.direction & INV))
            adjacent.visited = true;
        node1.edges.push_back(adjacent);
    }

    void ItstGraph::addNode(Node& node)
    {
        node.id = m_nodes.size();
        m_nodes.push_back(node);
        m_maps[node.vproxy.pointer()->idx] = node.id;
    }

    ItstGraph::ItstGraph(FH fh, ItstAlg* data, int meshId):
        mp_alg(data), m_fh(fh), m_meshId(meshId)
    {
        m_bValid = false;

        if (fh->isSimple())
        {
            for (int i = 0; i < 3; i++)
            {
                if (fh->edges[i]->data && fh->edges[i]->data->vertices.size())
                    m_bValid = true;
            }

            if (!m_bValid)
            {
                ReportError("Cannot build a graph!");
                return;
            }
        }

        m_bValid = true; 

        // 添加三个角点
        Node node;
        for (size_t i = 0; i < 3; i++)
        {
            if (!fh->vertices[i]->data)
            {
                assert(fh->data->planeRep);
                auto proxy = mp_alg->addVEntity(fh->data->planeRep->point(i));
                proxy.pointer()->addContext(meshId, fh, vertex_tag(i));
                fh->vertices[i]->data.reset(new UserVData(proxy));
            }
            node.vproxy = *fh->vertices[i]->data->proxy;
            addNode(node);
        }

        // 添加边界点, 以及相关边
        Edge edge;
        for (size_t i = 0; i < 3; i++)
        {
            if (!fh->edges[i]->data)
            {
                edge.vId[0] = fh->vertices[(i + 1) % 3]->data->proxy->pointer()->idx;
                edge.vId[1] = fh->vertices[(i + 2) % 3]->data->proxy->pointer()->idx;
                edge.direction = SEQ;
                addEdge(edge);
                continue;
            }

            std::deque<VProxyItr> points;
            for (size_t j= 0; j < fh->edges[i]->data->vertices.size(); j++)
            {
                points.push_back(fh->edges[i]->data->vertices[j]);
                node.vproxy = fh->edges[i]->data->vertices[j];
                addNode(node);
            }
            PBPointCompare cmp(fh->data->sp, fh->data->planeRep->m_bps[i]);
            std::sort(points.begin(), points.end(), cmp);

            edge.direction = SEQ;
            edge.vId[0] = fh->vertices[(i + 1) % 3]->data->proxy->pointer()->idx;
            edge.vId[1] = points[0].pointer()->idx;
            addEdge(edge);

            edge.vId[0] = points[points.size()-1].pointer()->idx;
            edge.vId[1] = fh->vertices[(i + 2) % 3]->data->proxy->pointer()->idx;
            addEdge(edge);

            for (size_t j = 0; j < points.size()-1; j++)
            {
                edge.vId[0] = points[j].pointer()->idx;
                edge.vId[1] = points[j+1].pointer()->idx;
                addEdge(edge);
            }
        }

        // 添加内点，注意这里没有考虑相交的情形
        if (!fh->data->itstTri || !fh->data->itstTri->isectLines.size()) return;

        auto &inners = fh->data->itstTri->inVertices;
        for (size_t i = 0; i < inners.size(); i++)
        {
            node.vproxy = inners[i];
            addNode(node);
        }

        auto &biedges = fh->data->itstTri->isectLines;
        edge.direction = SEQ | INV;
        for (auto &e: biedges)
        {
            edge.vId[0] = e.pts[0].gIdx;
            edge.vId[1] = e.pts[1].gIdx;
            addEdge(edge);
        }

        sortEdgeDirection();
    }


    ItstGraph::~ItstGraph()
    {
    }

    int determinePhase(double abcos, double absin)
    {
        if (abcos >= 0.0 && absin >= 0.0)
            return 1;

        if (abcos < 0.0 && absin >= 0.0)
            return 2;

        if (abcos < 0.0 && absin < 0.0)
            return 3;

        else return 4;
    }

    struct AngleComparison
    {
        typedef ItstGraph::Adj Item;

        AngleComparison(){}

        bool operator()(Item& a, Item&b)
        {
            auto v1 = a.calcDirection(graph, origin);
            auto v2 = b.calcDirection(graph, origin);

            double cos1 = v1 * zero;
            double cos2 = v2 * zero;

            double sin1 = sqrt(CGAL::cross_product(zero, v1).squared_length());
            double sin2 = sqrt(CGAL::cross_product(zero, v2).squared_length());

            int phase1 = determinePhase(cos1, sin1);
            int phase2 = determinePhase(cos2, sin2);

            if (phase2 < phase1) return false;
            if (phase2 > phase1) return true;

            double tan1 = sin1 / cos1;
            double tan2 = sin2 / cos2;

            return tan2 > tan1;
        }

        K::Vector_3 zero;
        PBPoint<K> origin;
        ItstGraph* graph = nullptr;
    };


    K::Vector_3 ItstGraph::Adj::calcDirection(ItstGraph* graph, const PBPoint<K>& origin)
    {
        auto there = graph->get_nodes()[vId].vproxy.pointer()->pos.getCoord();
        return there - origin.getCoord();
    }

    void ItstGraph::sortEdgeDirection()
    {
        for (auto &node : m_nodes)
        {
            // 以第一个为O点，计算sin, cos, 象限排序
            AngleComparison comp;
            comp.origin = node.vproxy.pointer()->pos;
            comp.zero = node.edges.begin()->calcDirection(this, comp.origin);
            comp.graph = this;
            std::sort(node.edges.begin() + 1, node.edges.end(), comp);
        }
    }

    void ItstGraph::floodFilling(int startId)
    {
        std::queue<int, std::list<int>> queue;
        m_nodes[startId].mark = SEEDED;
        queue.push(startId);

        while (!queue.empty())
        {
            startId = queue.front();
            auto &startNode = m_nodes[startId];

            if (startNode.mark == VISITED)
            {
                queue.pop();
                continue;
            }

            auto &adjs = startNode.edges;
            startNode.mark = VISITED;

            for (auto& edge : adjs)
            {
                if (isNodeVisited(edge)) continue;

                auto &node = m_nodes[edge.vId];
                node.indicator = new SampleIndicatorVector(m_ids);

                for (int id : m_ids)
                {
                    assert((*startNode.indicator)[id] != REL_UNKNOWN &&
                        (*startNode.indicator)[id] != REL_NOT_AVAILABLE);

                    if (!startNode.vproxy.pointer()->hasContext(id))
                        node.indicator->at(id) = startNode.indicator->at(id);
                }

                for (auto& ctx : node.vproxy.pointer()->ctx)
                    if (m_meshId != ctx.meshId)
                        node.indicator->at(ctx.meshId) = REL_ON_BOUNDARY;

                for (int id : m_ids)
                {
                    if (node.indicator->at(id) == REL_UNKNOWN)
                    {
                        auto location = startNode.vproxy.pointer()->findInContext(id);
                        node.indicator->at(id) = relationOfContext(*location, node.vproxy.pointer()->pos);
                    }
                }

                queue.push(edge.vId);
                node.mark = SEEDED;

            }
            queue.pop();
        }
    }

    void ItstGraph::getAllLoops(std::deque<Loop>& loops)
    {
        // 从第一个node开始，必然是一个角点
        typedef std::vector<Adj>::iterator AdjItr;
        std::queue<Adj*, std::list<Adj*>> queue;

        for (auto& adjitr : m_nodes[0].edges)
        {
            if (!adjitr.visited)
            {
                queue.push(&adjitr);
                break;
            }
        }

        while (!queue.empty())
        {
            if (!queue.front()->visited)
            {
                Loop loop;

                Adj* adj = queue.front();
                int startId = adj->startId;
                int begin = startId;
                int next = adj->vId;
                adj->visited = true;
                loop.push_back(&m_nodes[startId]);

                while (startId != next)
                {
                    auto& nextNode = m_nodes[next];
                    loop.push_back(&nextNode);

                    size_t n = nextNode.edges.size();
                    bool got = false;
                    for (size_t i = 0; i < n; i++)
                    {
                        if (nextNode.edges[i].vId == begin)
                        {
                            got = true;
                            adj = &nextNode.edges[(i + 1) % n];
                            begin = next;
                            next = adj->vId;
                            assert(!adj->visited);
                            adj->visited = true;

                            // add opposite if needed
                            if (!nextNode.edges[i].visited)
                                queue.push(&nextNode.edges[i]);

                            break;
                        }
                    }
                    assert(got);
                }
                loops.push_back(loop);
            }

            queue.pop();
        }
    }

}
