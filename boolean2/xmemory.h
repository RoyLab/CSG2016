#pragma once
#include <vector>
#include <macros.h>
#include "preps.h"
#include "global.h"

namespace Boolean
{
    typedef Depick CGALKernel;
    typedef typename CGALKernel::Triangle_3 CGALTriangle;
    typedef typename CGALKernel::Point_3 CGALPoint;

    struct MyVertex
    {
    public:
        typedef uint32_t Index;

    private:
        int rep = 0; // positive is v-base, negative is p-base
    public:
        void setAsPRep(int i) { rep = -(i + 1); }
        void setAsVRep(int i) { rep = (i + 1); }

        std::list<uint32_t> edges;

        bool findEdge(uint32_t other, uint32_t* result = nullptr) const;
        Oriented_side orientation(const XPlane& p) const;

        bool isPlaneRep() const { return rep < 0; }
        uint32_t id() const { assert(rep);  return std::abs(rep) - 1; }
        bool operator==(const XPoint& p) const;
        bool isValid() const { return rep != 0; }

        const XPoint& ppoint() const;
        const cyPointT& point() const;
    };

    struct FH
    {
        FH() {}
        FH(int o, IPolygon* p) : orientation(o), ptr(p) {}
        int orientation; // +1 same, -1 oppo
        IPolygon* ptr = nullptr;
    };

    struct MyEdge
    {
    public:
        typedef uint32_t Index;
        typedef int32_t SIndex;

        class ConstFaceIterator
        {
        public:
            ConstFaceIterator(const MyEdge& edge) :
                m_edge(edge), stage(0) {}

            ConstFaceIterator& operator++();
            operator bool() const { return stage != -1; }

            const FH& faceHandle() const
            {
                if (stage < 0) throw std::exception();
                if (stage < 2) return m_edge.fhs[stage];
                else return *eItr;
            }
            const IPolygon* face() const { return faceHandle().ptr; }
            int orientation() const { return faceHandle().orientation; }

        private:
            const MyEdge& m_edge;
            int stage; // > 2, then eItr
            std::list<FH>::const_iterator eItr;
        };

        class FaceIterator
        {
        public:
            FaceIterator(MyEdge& edge, bool triangle = false);

            FaceIterator& operator++();
            FaceIterator& incrementToTriangle();

            operator bool() const { return stage != -1; }
            FH& faceHandle() const
            {
                if (stage < 0) throw std::exception();
                if (stage < 2) return m_edge.fhs[stage];
                else return *eItr;
            }
            IPolygon* face() const { return faceHandle().ptr; }
            int orientation() const { return faceHandle().orientation; }

        private:
            MyEdge& m_edge;
            int stage; // > 2, then eItr
            std::list<FH>::iterator eItr;
        };

    public:
        MyVertex::Index ends[2];
        std::vector<NeighborInfo>* neighbor = nullptr;
        EdgeInsctData* inscts = nullptr;
        bool noOverlapNeighbor = false;

        MyEdge(MyVertex::Index a, MyVertex::Index b) : ends{ a, b } {}
        ~MyEdge();
        void addAjacentFace(MyVertex::Index s, MyVertex::Index e, IPolygon* fPtr);
        int faceOrientation(const IPolygon*) const;
        bool remove(const IPolygon*); // BUG: vertex √ª¥¶¿Ì
        uint32_t faceCount() const;
        MyVertex::Index theOtherVId(MyVertex::Index thiz) const;
        MyVertex& theOtherVertex(MyVertex::Index thiz) const;

    protected:
        FH fhs[2];
        std::list<FH> extrafhs;
    };

    class MemoryManager
    {
        typedef cyPointT VPoint;
    public:
        std::vector<XPlaneBase> planes;
        std::vector<MyEdge> edges;
        std::vector<MyVertex> vertices;

        std::vector<VPoint>	points;
        std::vector<XPoint>	ppoints;
        ExternPtr std::vector<Triangle*> insctTris;
        std::vector<SubPolygon*> subpolys;

    public:
        ~MemoryManager() {}
        static MemoryManager* getInstance();

        uint32_t insertVertices(VPoint* begin, VPoint* end);
        uint32_t insertVertex(XPoint& pt);
        uint32_t getEdgeId(uint32_t a, uint32_t b, IPolygon* facePtr);
        void addSubPolygon(SubPolygon* poly) { subpolys.push_back(poly); }
        void outputIntersection(const std::string&, const cyPointT&, const cyPointT&);
        void clear();

    private:
        MemoryManager() {}
        static MemoryManager mgr;
    };

    inline const MyEdge& xcedge(uint32_t id) { return MemoryManager::getInstance()->edges[id]; }
    inline MyEdge& xedge(uint32_t id) { return MemoryManager::getInstance()->edges[id]; }
    inline std::vector<MyEdge>& xedges() { return MemoryManager::getInstance()->edges; }

    inline const MyVertex& xcvertex(uint32_t id) { return MemoryManager::getInstance()->vertices[id]; }
    inline MyVertex& xvertex(uint32_t id) { return MemoryManager::getInstance()->vertices[id]; }
    inline std::vector<MyVertex>& xvertices() { return MemoryManager::getInstance()->vertices; }

    inline const cyPointT& xcpoint(uint32_t id) { return MemoryManager::getInstance()->points[id]; }
    inline cyPointT& xpoint(uint32_t id) { return MemoryManager::getInstance()->points[id]; }
    inline std::vector<cyPointT>& xpoints() { return MemoryManager::getInstance()->points; }

    inline const XPoint& xcppoint(uint32_t id) { return MemoryManager::getInstance()->ppoints[id]; }
    inline XPoint& xppoint(uint32_t id) { return MemoryManager::getInstance()->ppoints[id]; }
    inline std::vector<XPoint>& xppoints() { return MemoryManager::getInstance()->ppoints; }

    inline const XPlaneBase& xcplane(uint32_t id) { return MemoryManager::getInstance()->planes[id]; }
    inline XPlaneBase& xplane(uint32_t id) { return MemoryManager::getInstance()->planes[id]; }
    inline std::vector<XPlaneBase>& xplanes() { return MemoryManager::getInstance()->planes; }

    inline std::vector<Triangle*>& intersectTriangles() { return MemoryManager::getInstance()->insctTris; }

    // other functions
    int linearOrder(const XLine& l, const MyVertex& a, const MyVertex& b);
    static inline int linearOrder(const XLine& l, MyVertex::Index a, MyVertex::Index b) { return linearOrder(l, xvertex(a), xvertex(b)); }

    Oriented_side orientation(const XPlane& p, const MyVertex& v);
    static inline Oriented_side orientation(const XPlane& p, MyVertex::Index v) { return orientation(p, xvertex(v)); }
    XPlane pickPositiveVertical(const XLine&, const MyVertex&);
}
