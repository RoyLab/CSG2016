#pragma once
#include <vector>
#include <deque>
#include <set>

#include <macros.h>
#include "preps.h"
#include "global.h"
#include "RegularMesh.h"

namespace Boolean
{
    enum VertexType
    {
        V_Plane, V_PlanePoint, VI_Point
    };

    struct CsgOption
    {
        bool triangulation = false;
    };

    struct MergedVertex
    {
        std::set<VertexIndex> refs;
        std::vector<EdgeIndex> edges;
    };

    class GlobalData
    {
    public:
        std::deque<XPlaneBase>      planebases;

        std::deque<PlanePoint>	    ppoints; // iterator must be always valid
        std::deque<cyPointT>	    points; // iterator must be always valid

        std::deque<MyEdge>         edges;
        std::deque<MyVertex>       vertices;
        std::vector<RegularMesh*>   meshes;

        static GlobalData* getObject();

        int insertVertices(cyPointT* begin, cyPointT* end);
        VertexIndex insertVertex(PlanePoint& pt);
        EdgeIndex getEdgeId(VertexIndex a, VertexIndex b, IPolygon* facePtr);
        const std::vector<EdgeIndex>& get_merged_edges(VertexIndex id) const;
        void add_merged_edges(VertexIndex mergeid, EdgeIndex edgeid);
        void mergeVertices(VertexIndex a, VertexIndex b);

        void clear();

        // debug only
        void dumpIntersectionToXyzFile(const std::string& filename, const cyPointT& center, const cyPointT& scale);

    private:
        std::deque<MergedVertex> mergedvertices_;

        GlobalData() {}
        static GlobalData mgr;
    };

    /// useful functions
    void mergeVertices(std::set<VertexIndex>& indices);
    inline bool vertex_id_equals(VertexIndex a, VertexIndex b)
    {
        if (a == b) return true;
        return xvertex(a).id_equals(xvertex(b));
    }

    /// gramma sugar
    inline const MyEdge& xcedge(EdgeIndex id) { return GlobalData::getObject()->edges[id]; }
    inline MyEdge& xedge(EdgeIndex id) { return GlobalData::getObject()->edges[id]; }
    //inline std::vector<MyEdge>& xedges() { return GlobalData::getObject()->edges; }
    EdgeIndex assign_new_edge(MyEdge** p_edge);

    inline const MyVertex& xcvertex(VertexIndex id) { return GlobalData::getObject()->vertices[id]; }
    inline MyVertex& xvertex(VertexIndex id) { return GlobalData::getObject()->vertices[id]; }
    //inline std::vector<MyVertex>& xvertices() { return GlobalData::getObject()->vertices; }
    //inline VertexIndex assign_new_edge(MyVertex** p_edge);

    inline const cyPointT& xcpoint(Index id) { return GlobalData::getObject()->points[id]; }
    inline cyPointT& xpoint(Index id) { return GlobalData::getObject()->points[id]; }
    //inline std::deque<cyPointT>& xpoints() { return GlobalData::getObject()->points; }
    //inline EdgeIndex assign_new_edge(MyEdge** p_edge);

    inline const PlanePoint& xcppoint(Index id) { return GlobalData::getObject()->ppoints[id]; }
    inline PlanePoint& xppoint(Index id) { return GlobalData::getObject()->ppoints[id]; }
    //inline std::deque<PlanePoint>& xppoints() { return GlobalData::getObject()->ppoints; }
    //inline EdgeIndex assign_new_edge(MyEdge** p_edge);

    inline const XPlaneBase& xcplane(Index id) { return GlobalData::getObject()->planebases[id]; }
    inline XPlaneBase& xplane(Index id) { return GlobalData::getObject()->planebases[id]; }
    //inline std::deque<XPlaneBase>& xplanes() { return GlobalData::getObject()->planebases; }
    Index assign_new_plane(XPlaneBase** p_edge);

    inline const std::vector<RegularMesh*>& xcmeshlist() { return GlobalData::getObject()->meshes; }
    inline std::vector<RegularMesh*>& xmeshlist() { return GlobalData::getObject()->meshes; }

    /// more  gramma sugar
    //inline int linear_order(const PlaneLine& l, VertexIndex a, VertexIndex b) { return linear_order(l, xvertex(a), xvertex(b)); }
    //inline Oriented_side orientation(const XPlane& p, VertexIndex v) { return orientation(p, xvertex(v)); }

}
