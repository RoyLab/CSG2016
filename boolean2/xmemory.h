#pragma once
#include <vector>
#include <macros.h>
#include <deque>
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
        std::vector<PlanePoint>	    ppoints;

        std::vector<cyPointT>	    points;

        std::vector<MyEdge>         edges;
        std::vector<MyVertex>       vertices;
        std::vector<RegularMesh*>   meshes;

        static GlobalData* getObject();

        uint32_t insertVertices(cyPointT* begin, cyPointT* end);
        uint32_t insertVertex(PlanePoint& pt);
        uint32_t getEdgeId(uint32_t a, uint32_t b, IPolygon* facePtr);

        void clear();

        // debug only
        void dumpIntersectionToXyzFile(const std::string& filename, const cyPointT& center, const cyPointT& scale);

    private:
        GlobalData() {}
        static GlobalData mgr;
    };

    /// useful functions
    void mergeBrepVertices(VertexIndex a, VertexIndex b);

    /// gramma sugar
    inline const MyEdge& xcedge(uint32_t id) { return GlobalData::getObject()->edges[id]; }
    inline MyEdge& xedge(uint32_t id) { return GlobalData::getObject()->edges[id]; }
    inline std::vector<MyEdge>& xedges() { return GlobalData::getObject()->edges; }

    inline const MyVertex& xcvertex(uint32_t id) { return GlobalData::getObject()->vertices[id]; }
    inline MyVertex& xvertex(uint32_t id) { return GlobalData::getObject()->vertices[id]; }
    inline std::vector<MyVertex>& xvertices() { return GlobalData::getObject()->vertices; }

    inline const cyPointT& xcpoint(uint32_t id) { return GlobalData::getObject()->points[id]; }
    inline cyPointT& xpoint(uint32_t id) { return GlobalData::getObject()->points[id]; }
    inline std::vector<cyPointT>& xpoints() { return GlobalData::getObject()->points; }

    inline const PlanePoint& xcppoint(uint32_t id) { return GlobalData::getObject()->ppoints[id]; }
    inline PlanePoint& xppoint(uint32_t id) { return GlobalData::getObject()->ppoints[id]; }
    inline std::vector<PlanePoint>& xppoints() { return GlobalData::getObject()->ppoints; }

    inline const XPlaneBase& xcplane(uint32_t id) { return GlobalData::getObject()->planebases[id]; }
    inline XPlaneBase& xplane(uint32_t id) { return GlobalData::getObject()->planebases[id]; }
    inline std::deque<XPlaneBase>& xplanes() { return GlobalData::getObject()->planebases; }

    inline const std::vector<RegularMesh*>& xcmeshlist() { return GlobalData::getObject()->meshes; }
    inline std::vector<RegularMesh*>& xmeshlist() { return GlobalData::getObject()->meshes; }

    /// more  gramma sugar
    //inline int linear_order(const PlaneLine& l, VertexIndex a, VertexIndex b) { return linear_order(l, xvertex(a), xvertex(b)); }
    //inline Oriented_side orientation(const XPlane& p, VertexIndex v) { return orientation(p, xvertex(v)); }

}
