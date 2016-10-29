#include "precompile.h"
#include "xmemory.h"
#include "RegularMesh.h"


namespace Boolean
{
    MemoryManager * MemoryManager::getInstance()
    {
        static MemoryManager mgr;
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
        return ppoints.size() - 1;
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
}
