#include "precompile.h"
#include <xgeometry.h>

#include "xmemory.h"
#include "RegularMesh.h"


namespace Boolean
{
    GlobalData GlobalData::mgr;

    GlobalData * GlobalData::getObject()
    {
        return &mgr;
    }

    uint32_t GlobalData::insertVertices(cyPointT * begin, cyPointT * end)
    {
        int offset = points.size();
        int n = end - begin;

        points.insert(points.end(), begin, end);
        vertices.resize(vertices.size() + n);

        for (int i = 0; i < n; i++)
            vertices[offset + i].setAsVRep(offset + i);

        return offset;
    }

    uint32_t GlobalData::insertVertex(PlanePoint & pt)
    {
        ppoints.push_back(pt);
        int vId =  ppoints.size() - 1;
        MyVertex ver;
        ver.setAsPRep(vId);
        vertices.push_back(ver);
        return vertices.size() - 1;
    }

    uint32_t GlobalData::getEdgeId(uint32_t a, uint32_t b, IPolygon * facePtr)
    {
        assert(a < xvertices().size());
        MyVertex& one = xvertex(a);

        uint32_t target;
        bool res = one.findEdge(b, &target);

        if (!res)
        {
            MyVertex& theother = xvertex(b);
            assert(!theother.findEdge(a));
            
            target = edges.size();
            edges.emplace_back(a, b);

            one.edges.push_back(target);
            theother.edges.push_back(target);
        }
        xedge(target).addAjacentFace(a, b, facePtr);
        return target;
    }

    void GlobalData::dumpIntersectionToXyzFile(const std::string &fileName, const cyPointT& center, const cyPointT& scale)
    {
        std::ofstream file(fileName);

        for (auto& p : ppoints)
        {
            auto p2 = p.toVertexBased();
            XR::invCoords(center, scale, p2);
            file << p2.x << ' ' << p2.y << ' ' << p2.z << std::endl;
        }

        file.close();
    }

    void GlobalData::clear()
    {
        planebases.clear();
        edges.clear();
        vertices.clear();
        points.clear();
        ppoints.clear();
    }

    MyEdge::ConstFaceIterator & MyEdge::ConstFaceIterator::operator++()
    {
        if (stage == -1) throw std::exception();

        if (stage < 2)
        {
            ++stage;
            if (stage == 2)
            {
                if (m_edge.extrafhs.empty())
                {
                    stage = -1;
                    return *this;
                }
                else eItr = m_edge.extrafhs.begin();
            }
        }
        else
        {
            ++eItr;
            if (eItr == m_edge.extrafhs.end())
                stage = -1;
        }
        return *this;
    }

    MyEdge::FaceIterator::FaceIterator(MyEdge & edge, bool triangle):
        m_edge(edge), stage(0)
    {
        if (triangle)
        {
            if (face()->getType() != IPolygon::TRIANGLE)
                incrementToTriangle();
        }
    }

    MyEdge::FaceIterator & MyEdge::FaceIterator::operator++()
    {
        if (stage == -1) throw std::exception();

        if (stage < 2)
        {
            ++stage;
            if (stage == 2)
            {
                if (m_edge.extrafhs.empty())
                {
                    stage = -1;
                    return *this;
                }
                else eItr = m_edge.extrafhs.begin();
            }
        }
        else
        {
            ++eItr;
            if (eItr == m_edge.extrafhs.end())
                stage = -1;
        }
        return *this;
    }

    MyEdge::FaceIterator & MyEdge::FaceIterator::incrementToTriangle()
    {
        do {
        ++*this;
        } while (*this && face()->getType() != IPolygon::TRIANGLE);
        return *this;
    }

    void mergeBrepVertices(VertexIndex to, VertexIndex from)
    {
        MyEdge& edge1 = fh[to]->edge((from + 1) % 3);
        MyEdge& edge2 = fh[to]->edge((from + 2) % 3);

        if (edge1.ends[0] == *slots[to])
        {
            edge1.ends[0] = vid;
        }
        else
        {
            if (edge1.ends[1] == *slots[to])
                edge1.ends[1] = vid;
        }

        if (edge2.ends[0] == *slots[to])
        {
            edge2.ends[0] = vid;
        }
        else
        {
            if (edge2.ends[1] == *slots[to])
                edge2.ends[1] = vid;
        }

        MyVertex& vRef = xvertex(*slots[(to + 1) % 2]),
            &vMerge = xvertex(*slots[to]);

        if (!vMerge.edges().empty())
        {
            vRef.edges().insert(vRef.edges().end(), vMerge.edges().begin(), vMerge.edges().end());
            vMerge.edges().clear();
        }
    }
}
