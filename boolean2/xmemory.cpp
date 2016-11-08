#include "precompile.h"
#include "xmemory.h"
#include "RegularMesh.h"


namespace Boolean
{
    MemoryManager MemoryManager::mgr;

    MemoryManager * MemoryManager::getInstance()
    {
        return &mgr;
    }

    uint32_t MemoryManager::insertVertices(VPoint * begin, VPoint * end)
    {
        int offset = points.size();
        int n = end - begin;

        points.insert(points.end(), begin, end);
        vertices.resize(vertices.size() + n);

        for (int i = 0; i < n; i++)
            vertices[offset + i].setAsVRep(offset + i);

        return offset;
    }

    uint32_t MemoryManager::insertVertex(XPoint & pt)
    {
        ppoints.push_back(pt);
        int vId =  ppoints.size() - 1;
        MyVertex ver;
        ver.setAsPRep(vId);
        vertices.push_back(ver);
        return vertices.size() - 1;
    }

    uint32_t MemoryManager::getEdgeId(uint32_t a, uint32_t b, IPolygon * facePtr)
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

    void MemoryManager::outputIntersection(const std::string &fileName, const cyPointT& center, const cyPointT& scale)
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

}
