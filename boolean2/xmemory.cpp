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

    int GlobalData::insertVertices(cyPointT * begin, cyPointT * end)
    {
        int offset = points.size();
        int n = end - begin;

        points.insert(points.end(), begin, end);
        vertices.resize(vertices.size() + n);

        for (int i = 0; i < n; i++)
            vertices[offset + i].setAsVRep(offset + i);

        return offset;
    }

    VertexIndex GlobalData::insertVertex(PlanePoint & pt)
    {
        ppoints.push_back(pt);
        int vId =  ppoints.size() - 1;
        MyVertex ver;
        ver.setAsPRep(vId);
        vertices.push_back(ver);
        return vertices.size() - 1;
    }

    EdgeIndex GlobalData::getEdgeId(VertexIndex a, VertexIndex b, IPolygon * facePtr)
    {
        assert(a < vertices.size());
        MyVertex& one = xvertex(a);

        EdgeIndex target;
        bool res = one.findEdge(b, &target);

        if (!res)
        {
            MyVertex& theother = xvertex(b);
            assert(!theother.findEdge(a));
            
            target = edges.size();
            edges.emplace_back(a, b);

            one.add_edge(target);
            theother.add_edge(target);
        }
        xedge(target).addAjacentFace(a, b, facePtr);
        return target;
    }

    const std::vector<EdgeIndex>& GlobalData::get_merged_edges(VertexIndex id) const
    {
        return mergedvertices_[id].edges;
    }

    void GlobalData::add_merged_edges(VertexIndex mergeid, EdgeIndex edgeid)
    {
        mergedvertices_[mergeid].edges.push_back(edgeid);
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

    //void mergeVertices(VertexIndex to, VertexIndex from)
    //{
        //MyEdge& edge1 = fh[to]->edge((from + 1) % 3);
        //MyEdge& edge2 = fh[to]->edge((from + 2) % 3);

        //if (edge1.ends[0] == *slots[to])
        //{
        //    edge1.ends[0] = vid;
        //}
        //else
        //{
        //    if (edge1.ends[1]   == *slots[to])
        //        edge1.ends[1] = vid;
        //}

        //if (edge2.ends[0] == *slots[to])
        //{
        //    edge2.ends[0] = vid;
        //}
        //else
        //{
        //    if (edge2.ends[1] == *slots[to])
        //        edge2.ends[1] = vid;
        //}

        //MyVertex& vRef = xvertex(*slots[(to + 1) % 2]),
        //    &vMerge = xvertex(*slots[to]);

    void mergeVertices(VertexIndex a, VertexIndex b)
    {

    }

    void mergeVertices(std::set<VertexIndex>& indices)
    {
        VertexIndex v0 = INVALID_UINT32;
        for (VertexIndex v: indices)
        {
            if (v0 = INVALID_UINT32)
            {
                v0 = v;
            }
            else
            {
                mergeVertices(v0, v);
            }
        }
    }

    Index assign_new_plane(XPlaneBase ** p_edge)
    {
        auto &container = GlobalData::getObject()->planebases;
        container.emplace_back();
        *p_edge = &container.back();
        return container.size() - 1;
    }

    //if (!vMerge.edges().empty())
        //{
        //    vRef.edges().insert(vRef.edges().end(), vMerge.edges().begin(), vMerge.edges().end());
        //    vMerge.edges().clear();
        //}
    //}
}
