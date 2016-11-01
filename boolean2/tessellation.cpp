#include "precompile.h"
#include <list>
#include <map>
#include "xmemory.h"
#include "RegularMesh.h"
#include "intersection.h"
#include "XStruct.hpp"
#include "XException.hpp"

namespace Boolean
{
    struct PlaneVertex
    {
        XPlane plane;
        MyVertex::Index id;
    };

    struct EdgeAuxiliaryStructure
    {
        MyVertex::Index start, end;
        XLine line;
    };

    int linearOrder(const XLine& l, const MyVertex& a, const MyVertex& b)
    {
        int type = 0;
        if (a.isPlaneRep()) type += 1;
        if (b.isPlaneRep()) type += 2;

        switch (type)
        {
        case 0:
            return l.linearOrder(a.point(), b.point()) > 0;
        case 1:
            return l.pickPositiveVertical(a.ppoint())
                .orientation(b.point()) == ON_POSITIVE_SIDE;
        case 2:
            return l.pickPositiveVertical(b.ppoint())
                .orientation(a.point()) == ON_NEGATIVE_SIDE;
        case 3:
            return l.linearOrder(a.ppoint(), b.ppoint()) > 0;
        default:
            throw std::exception();
        }
        return 0;
    }


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

    void InsctData<EdgePBI>::refine(void* pData)
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
        for (int i = 0; i < seqs.size(); i++)
        {
            idmap[seqs[i].id] = i;
            newPbi[i].ends[1] = seqs[i].id;
            newPbi[i+1].ends[0] = seqs[i].id;
            newPbi[i].pends[1] = seqs[i].plane;
            newPbi[i + 1].pends[0] = seqs[i].plane;
        }
        newPbi[0].ends[0] = data.start;
        newPbi[points.size()].ends[1] = data.end;

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

    void InsctData<FacePBI>::refine(void* pData)
    {
        if (isRefined() || inscts.size() < 2) return;

        bRefined = true;
        throw XR::NotImplementedException();
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

                if (!e.inscts)
                {
                    //assert(m_nodes.find(e.ends[0]) != m_nodes.end());
                    //assert(m_nodes.find(e.ends[1]) != m_nodes.end());
                    //edge.v[0] = m_nodes.find(e.ends[0]);
                    //edge.v[1] = m_nodes.find(e.ends[1]);
                    //m_edges.push_back(edge);

                    continue;
                }

                int ires = e.faceOrientation(tri);
                assert(ires != 0);
                int i0 = 0, i1 = 1;
                if (ires < 0) std::swap(i0, i1);

                for (auto &set : e.inscts->inscts)
                {
                    for (auto& ePBI : set.second)
                    {
                        assert(m_nodes.find(ePBI.ends[0]) != m_nodes.end());
                        assert(m_nodes.find(ePBI.ends[1]) != m_nodes.end());

                        edge.v[i0] = m_nodes.find(ePBI.ends[0]);
                        edge.v[i1] = m_nodes.find(ePBI.ends[1]);
                        assert(edge.checkPlaneOrientation(tri));

                        m_edges.push_back(edge);

                        assert(edge.v[0] != m_nodes.end());
                        assert(edge.v[1] != m_nodes.end());

                        edge.v[i0]->second.edges.push_back(m_edges.size() - 1);
                        edge.v[i1]->second.edges.push_back(m_edges.size() - 1);
                    }
                }
            }

            edge.dir = D_NODIR;
            for (auto &set : tri->inscts->inscts)
            {
                for (auto& fPBI : set.second)
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
        bool TessGraph::Edge::checkPlaneOrientation(const Triangle *pTri)
        {
            XLine line(pTri->supportingPlane(), prep);
            if (linearOrder(line, xvertex(v[0]->first),
                xvertex(v[1]->first)) > 0)
                return true;

            return false;
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
                pTri->inscts->refine((void*)&pTri->supportingPlane());

            for (int i = 0; i < 3; i++)
            {
                EdgeAuxiliaryStructure data = { pTri->edge(i).ends[0] , pTri->edge(i).ends[1] };
                if (pTri->edge(i).ends[0] != pTri->vertexId((i+1)%3))
                    data.line = XLine(pTri->supportingPlane(), pTri->boundingPlane(i).opposite());
                else
                    data.line = XLine(pTri->supportingPlane(), pTri->boundingPlane(i));

                 if (pTri->edge(i).inscts)
                    pTri->edge(i).inscts->refine((void*)&data);
            }

            TessGraph tg(pTri);
            tg.tessellate();
        }
    }
}