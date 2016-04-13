#include <CGAL\Cartesian_matrix.h>
#include <queue>
#include <list>
#include "ItstGraph.h"


namespace CSG
{
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

            auto cross1 = CGAL::cross_product(zero, v1);
            auto cross2 = CGAL::cross_product(zero, v2);

            double sign1 = cross1 * normal;
            double sign2 = cross2 * normal;

            int phase1 = determinePhase(cos1, sign1);
            int phase2 = determinePhase(cos2, sign2);

            if (phase1 > phase2) return true;
            if (phase2 > phase1) return false;

            double tan1 = sign1 / cos1;
            double tan2 = sign2 / cos2;

            return tan1 > tan2;
        }

        K::Vector_3 zero, normal;
        PBPoint<K> origin;
        ItstGraph* graph = nullptr;
    };


    K::Vector_3 ItstGraph::Adj::calcDirection(ItstGraph* graph, const PBPoint<K>& origin)
    {
        auto there = graph->get_nodes()[vId].vproxy.pointer()->pos.getCoord();
        return there - origin.getCoord();
    }

    ItstGraph::ItstGraph(FH fh, ItstAlg* data, int meshId, std::vector<MyMesh*>* meshes) :
        mp_alg(data), m_fh(fh), m_meshId(meshId), mp_meshList(meshes)
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
        if (!fh->data->itstTri || fh->data->itstTri->isectLines.empty()) return;

        auto tri = fh->data->itstTri;
        mp_alg->resolveIntersection(fh, meshId);

        auto &inners = tri->inVertices;
        for (size_t i = 0; i < inners.size(); i++)
        {
            node.vproxy = inners[i];
            addNode(node);
        }

        auto &biedges = fh->data->itstTri->unifiedLines;
        edge.direction = SEQ | INV;
        for (auto &e : biedges)
        {
            edge.vId[0] = e.pts[0].vertex.pointer()->idx;
            edge.vId[1] = e.pts[1].vertex.pointer()->idx;
            addEdge(edge);
        }

        sortEdgeDirection();
    }


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

    void ItstGraph::sortEdgeDirection()
    {
        for (auto &node : m_nodes)
        {
            // 以第一个为O点，计算sin, cos, 象限排序
            AngleComparison comp;
            comp.origin = node.vproxy.pointer()->pos;
            comp.zero = node.edges.begin()->calcDirection(this, comp.origin);
            comp.graph = this;
            comp.normal = m_fh->data->sp.orthogonal_vector();
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

                /* 如果种子没有context，则直接继承种子的这部分indicators */
                for (int id : m_ids)
                {
                    assert((*startNode.indicator)[id] != REL_UNKNOWN &&
                        (*startNode.indicator)[id] != REL_NOT_AVAILABLE);

                    if (!startNode.vproxy.pointer()->hasContext(id))
                        node.indicator->at(id) = startNode.indicator->at(id);
                }

                /* 如果衍生有context，直接为REL_ON_BOUNDARY */
                for (auto& ctx : node.vproxy.pointer()->ctx)
                    if (m_meshId != ctx.meshId)
                        node.indicator->at(ctx.meshId) = REL_ON_BOUNDARY;

                /* 剩下的种子有context而衍生没有的，用relationOfContext求解 */
                for (int id : m_ids)
                {
                    if (node.indicator->at(id) == REL_UNKNOWN)
                    {
                        auto location = startNode.vproxy.pointer()->findInContext(id);
                        auto relation = relationOfContext(*location, node.vproxy.pointer()->pos);

                        /* 修正结果 */
                        if (mp_meshList->at(id)->bInverse)
                            inverseRelation(relation);

                        node.indicator->at(id) = relation;
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

#ifdef _DEBUG
        for (auto& node : m_nodes)
        {
            for (auto& adj : node.edges)
                assert(adj.visited);
        }
#endif
    }

}
