#include "precompile.h"
#include <list>
#include <map>
#include "xmemory.h"
#include "RegularMesh.h"
#include "intersection.h"
#include "XStruct.hpp"

namespace Boolean
{
    void InsctData<EdgePBI>::refine()
    {
        if (isRefined()) return;

        bRefined = true;
    }

    void InsctData<FacePBI>::refine()
    {
        if (isRefined()) return;

        bRefined = true;
    }

    namespace
    {
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

                VertexIndex startVertex() const { return v[0]->first; }
                VertexIndex endVertex() const { return v[1]->first; }

                bool checkPlaneOrientation(const Triangle*);
            };

            struct SortObject
            {
                bool operator() (EdgeIndex i, EdgeIndex j);
                const TessGraph* pTG;
                VertexIndex node;
            };

        public:
            TessGraph(const Triangle*);
            void tessellate();

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
                for (auto vId : tri->inscts->points)
                    m_nodes[vId] = node;
            }

            edge.dir = D_SEQ;
            for (int i = 0; i < 3; i++)
            {
                auto& e = tri->edge(i);
                edge.prep = tri->boundingPlane(i);
                assert(edge.checkPlaneOrientation(tri));

                if (!e.inscts) continue;

                int ires = e.faceOrientation(tri);
                assert(ires != 0);
                int i0 = 0, i1 = 1;
                if (ires < 0) std::swap(i0, i1);

                for (auto &ePBI : e.inscts->inscts)
                {
                    assert(m_nodes.find(ePBI.ends[0]) != m_nodes.end());
                    assert(m_nodes.find(ePBI.ends[1]) != m_nodes.end());

                    edge.v[i0] = m_nodes.find(ePBI.ends[0]);
                    edge.v[i1] = m_nodes.find(ePBI.ends[1]);
                    m_edges.push_back(edge);

                    assert(edge.v[0] != m_nodes.end());
                    assert(edge.v[1] != m_nodes.end());

                    edge.v[i0]->second.edges.push_back(m_edges.size() - 1);
                    edge.v[i1]->second.edges.push_back(m_edges.size() - 1);
                }
            }

            edge.dir = D_NODIR;
            for (auto &fPBI : tri->inscts->inscts)
            {
                assert(m_nodes.find(fPBI.ends[0]) != m_nodes.end());
                assert(m_nodes.find(fPBI.ends[1]) != m_nodes.end());

                edge.v[0] = m_nodes.find(fPBI.ends[0]);
                edge.v[1] = m_nodes.find(fPBI.ends[1]);
                edge.prep = fPBI.vertPlane;
                assert(edge.checkPlaneOrientation(tri));

                m_edges.push_back(edge);
                edge.v[0]->second.edges.push_back(m_edges.size() - 1);
                edge.v[1]->second.edges.push_back(m_edges.size() - 1);
            }
        }

        void TessGraph::tessellate()
        {
            // sort all node by circular order
            SortObject sortObj = { this };
            for (auto& nPair : m_nodes)
            {
                auto &node = nPair.second;
                sortObj.node = nPair.first;
                std::quicksort(node.edges.begin(), node.edges.end(), sortObj);
                assert([&](std::vector<uint32_t>& vec)->bool {
                    for (int i = 1; i < vec.size(); i++)
                    {
                        if (!sortObj(vec[i - 1], vec[i]))
                            return false;
                    }
                    return true;
                }(node.edges));
            }

            // remove the original 
            bool bRes;
            for (int i = 0; i < 3; i++)
                bRes = xedge(mp_tri->edgeId(i)).remove(mp_tri);

            // tessellate
            typedef uint32_t VeretxIndex;
            NodeMap::iterator chead, cend;
            EdgeIndex cedge;
            std::stack<EdgeIndex> edgeStack;
            std::vector<uint32_t> loop;
            edgeStack.push(0);
            cend = m_edges[0].v[1];
            chead = m_edges[0].v[0];
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
                loop.push_back((eRef.dir == D_SEQ) ? eRef.v[1]->first : eRef.v[0]->first);

                cend = (eRef.dir == D_SEQ) ? eRef.v[1] : eRef.v[0];
                while (cend != chead)
                {
                    EdgeIndex newEId = cend->second.findnext(cedge);
                    Edge& newEdge = m_edges[newEId];

                    Direction dir = D_SEQ;
                    if (newEdge.v[0] != cend)
                        dir = D_INV;

                    assert(newEdge.dir & dir);
                    newEdge.dir = Direction(newEdge.dir - dir);
                    loop.push_back(cend->first);

                    if (newEdge.dir != D_NAN)
                        edgeStack.push(newEId);
                }

                // add subpolygon
                const uint32_t n = loop.size();
                SubPolygon *spoly = new SubPolygon(n);
                spoly->constructFromVertexList(loop.begin(), loop.end());
            }
        }

        bool TessGraph::SortObject::operator()(EdgeIndex i, EdgeIndex j)
        {
            auto &ei = pTG->m_edges[i];
            auto &ej = pTG->m_edges[j];

            auto sp = pTG->mp_tri->supportingPlane();
            XPlane ep[2] = { ei.prep, ej.prep };
            
            if (ei.startVertex() != node)
                ep[0].inverse();

            if (ej.startVertex() != node)
                ep[1].inverse();

            Real res = sign(sp, ep[0], ep[1]);
            assert(res != Real(0));
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
}

    void tessellation(std::vector<RegularMesh*>& meshes)
    {
        auto pMem = MemoryManager::getInstance();
        auto& inscts = intersectTriangles();
        for (Triangle* pTri : inscts)
        {
            assert(pTri->isAdded4Tess());
            if (pTri->inscts)
                pTri->inscts->refine();

            for (int i = 0; i < 3; i++)
                if (pTri->edge(i).inscts)
                    pTri->edge(i).inscts->refine();

            TessGraph tg(pTri);
            tg.tessellate();
        }
    }
}