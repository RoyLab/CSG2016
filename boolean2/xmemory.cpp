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

    EdgeIndex GlobalData::get_edge_id_local_scope(VertexIndex a, VertexIndex b, IPolygon * facePtr)
    {
        assert(a < vertices.size());
        MyVertex& one = xvertex(a);

        EdgeIndex target;
        bool res = one.findEdge_local(b, &target);

        if (!res)
        {
            MyVertex& theother = xvertex(b);
            assert(!theother.findEdge_local(a));
            
            target = edges.size();
            edges.emplace_back(a, b);

            one.add_edge_local(target);
            theother.add_edge_local(target);
        }
        xedge(target).addAjacentFace(a, b, facePtr);
        return target;
    }

    //const std::vector<EdgeIndex>& GlobalData::get_merged_edges(VertexIndex id) const
    //{
    //    return mergedvertices_[id].edges;
    //}

    //void GlobalData::add_merged_edges(VertexIndex mergeid, EdgeIndex edgeid)
    //{
    //    //mergedvertices_[mergeid].edges.push_back(edgeid);
    //}

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
        meshes.clear();
        mergedvertices_.clear();
        adj_graph.release();
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

    VertexIndex GlobalData::get_main_vertexId(VertexIndex v) const
    {
        if (xvertex(v).merge_ >= 0)
        {
            return *mergedvertices_[xvertex(v).merge_].refs.begin();
        }
        else
        {
            return v;
        }
    }

    void GlobalData::mergeVertices(VertexIndex a, VertexIndex b)
    {
        MyVertex& va = xvertex(a);
        MyVertex& vb = xvertex(b);

        if (va.merge_ < 0)
        {
            if (vb.merge_ < 0)
            {
                va.merge_ = vb.merge_ = mergedvertices_.size();

                mergedvertices_.emplace_back();
                MergedVertex& last = mergedvertices_.back();

                last.refs.insert(a);
                last.refs.insert(b);

                //last.edges.insert(last.edges.end(), va.edges_.begin(), va.edges_.end());
                //last.edges.insert(last.edges.end(), vb.edges_.begin(), vb.edges_.end());
            }
            else
            {
                va.merge_ = vb.merge_;
                MergedVertex& now = mergedvertices_[vb.merge_];

                now.refs.insert(a);
                //now.edges.insert(now.edges.end(), va.edges_.begin(), va.edges_.end());
            }
        }
        else
        {
            if (vb.merge_ < 0)
            {
                vb.merge_ = va.merge_;
                MergedVertex& now = mergedvertices_[va.merge_];

                now.refs.insert(b);
                //now.edges.insert(now.edges.end(), vb.edges_.begin(), vb.edges_.end());
            }
            else
            {
                if (va.merge_ == vb.merge_)
                {
                    return;
                }

                MergedVertex& now = mergedvertices_[va.merge_];
                MergedVertex& gotovanish = mergedvertices_[vb.merge_];

                now.refs.insert(gotovanish.refs.begin(), gotovanish.refs.end());
                //now.edges.insert(now.edges.end(), gotovanish.edges.begin(), gotovanish.edges.end());

                for (VertexIndex vidx : gotovanish.refs)
                {
                    xvertex(vidx).merge_ = va.merge_;
                }
            }
        }
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
                GlobalData::getObject()->mergeVertices(v0, v);
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
