#include <CGAL\Cartesian_matrix.h>

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
            Plane pa = a.planes[2];
            Plane pb = b.planes[2]; //之所以是2,因为2一定是个非支撑平面
            
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

            CGAL::Point_3<K> posa(aa.coord[0], aa.coord[1], aa.coord[2]);
            CGAL::Point_3<K> posb(bb.coord[0], bb.coord[1], bb.coord[2]);

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
        if (e.direction & SEQ)
        {
            adjacent.vId = node1.vproxy.pointer()->idx;
            node0.edges.push_back(adjacent);
        }

        if (e.direction & INV)
        {
            adjacent.vId = node0.vproxy.pointer()->idx;
            node1.edges.push_back(adjacent);
        }
    }

    void ItstGraph::addNode(Node& node)
    {
        m_nodes.push_back(node);
        m_maps[node.vproxy.pointer()->idx] = m_nodes.size() - 1;
    }

    ItstGraph::ItstGraph(FH fh, ItstAlg* data, int meshId):
        mp_alg(data)
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

        // 添加内点
        if (!fh->data->itstTri || !fh->data->itstTri->inVertices.size()) return;

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

    void ItstGraph::sortEdgeDirection()
    {

    }

}
