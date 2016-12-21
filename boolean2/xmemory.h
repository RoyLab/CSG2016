#pragma once
#include <vector>
#include <macros.h>
#include <deque>
#include <set>

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

    class GlobalData
    {
    public:
        std::deque<XPlaneBase>      planebases;

        std::deque<PlanePoint>	    ppoints; // iterator must be always valid
        std::deque<cyPointT>	    points; // iterator must be always valid

        std::vector<MyEdge>         edges;
        std::vector<MyVertex>       vertices;
        std::vector<RegularMesh*>   meshes;

        static GlobalData* getObject();

        int insertVertices(cyPointT* begin, cyPointT* end);
        VertexIndex insertVertex(PlanePoint& pt);
        EdgeIndex getEdgeId(VertexIndex a, VertexIndex b, IPolygon* facePtr);

        void clear();

        // debug only
        void dumpIntersectionToXyzFile(const std::string& filename, const cyPointT& center, const cyPointT& scale);

    private:
        GlobalData() {}
        static GlobalData mgr;
    };

    /// useful functions
    void mergeVertices(VertexIndex a, VertexIndex b);
    void mergeVertices(std::set<VertexIndex>& indices);

    /// gramma sugar
    inline const MyEdge& xcedge(EdgeIndex id) { return GlobalData::getObject()->edges[id]; }
    inline MyEdge& xedge(EdgeIndex id) { return GlobalData::getObject()->edges[id]; }
    inline std::vector<MyEdge>& xedges() { return GlobalData::getObject()->edges; }

    inline const MyVertex& xcvertex(VertexIndex id) { return GlobalData::getObject()->vertices[id]; }
    inline MyVertex& xvertex(VertexIndex id) { return GlobalData::getObject()->vertices[id]; }
    inline std::vector<MyVertex>& xvertices() { return GlobalData::getObject()->vertices; }

    inline const cyPointT& xcpoint(Index id) { return GlobalData::getObject()->points[id]; }
    inline cyPointT& xpoint(Index id) { return GlobalData::getObject()->points[id]; }
    inline std::deque<cyPointT>& xpoints() { return GlobalData::getObject()->points; }

    inline const PlanePoint& xcppoint(Index id) { return GlobalData::getObject()->ppoints[id]; }
    inline PlanePoint& xppoint(Index id) { return GlobalData::getObject()->ppoints[id]; }
    inline std::deque<PlanePoint>& xppoints() { return GlobalData::getObject()->ppoints; }

    inline const XPlaneBase& xcplane(Index id) { return GlobalData::getObject()->planebases[id]; }
    inline XPlaneBase& xplane(Index id) { return GlobalData::getObject()->planebases[id]; }
    inline std::deque<XPlaneBase>& xplanes() { return GlobalData::getObject()->planebases; }

    inline const std::vector<RegularMesh*>& xcmeshlist() { return GlobalData::getObject()->meshes; }
    inline std::vector<RegularMesh*>& xmeshlist() { return GlobalData::getObject()->meshes; }

    /// more  gramma sugar
    //inline int linear_order(const PlaneLine& l, VertexIndex a, VertexIndex b) { return linear_order(l, xvertex(a), xvertex(b)); }
    //inline Oriented_side orientation(const XPlane& p, VertexIndex v) { return orientation(p, xvertex(v)); }

}
