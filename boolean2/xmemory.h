#pragma once
#include <vector>
#include "preps.h"

namespace Boolean
{
    typedef Depick CGALKernel;
    typedef typename CGALKernel::Triangle_3 CGALTriangle;
    typedef typename CGALKernel::Point_3 CGALPoint;

    struct MyVertex
    {
    private:
        int rep = 0; // positive is v-base, negative is p-base
    public:
        void setAsPRep(int i) { rep = -(i + 1); }
        void setAsVRep(int i) { rep = (i + 1); }

        std::list<uint32_t> edges;

        bool findEdge(uint32_t other, uint32_t* result = nullptr) const;
        Oriented_side orientation(const XPlane&) const;

        bool isPlaneRep() const { return rep < 0; }
        uint32_t id() const { assert(rep);  return std::abs(rep) - 1; }
        bool operator==(const XPoint& p) const;
    };

    struct FH
    {
        FH() {}
        FH(int o, IPolygon* p) : orientation(o), ptr(p) {}
        int orientation;
        IPolygon* ptr = nullptr;
    };

    struct MyEdge
    {
    public:
        class ConstFaceIterator
        {
        public:
            ConstFaceIterator(const MyEdge& edge) :
                m_edge(edge), stage(0) {}

            ConstFaceIterator& operator++();
            ConstFaceIterator& operator++(int);
            operator bool() const;

            const IPolygon* ptr() const;

        private:
            const MyEdge& m_edge;
            int stage;
            std::list<FH>::const_iterator eItr;
        };
    public:
        MyEdge(uint32_t a, uint32_t b) : ends{ a, b } {}

        uint32_t ends[2];

        FH fhs[2];
        std::list<FH> extrafhs;
        InsctData<EdgePBI>* inscts = nullptr;

        ~MyEdge();
        void addAjacentFace(uint32_t s, uint32_t e, IPolygon* fPtr);
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

	public:
		~MemoryManager() {}
		static MemoryManager* getInstance();

        uint32_t insertVertices(VPoint* begin, VPoint* end);
        uint32_t insertVertex(XPoint& pt);
		uint32_t getEdgeId(uint32_t a, uint32_t b, IPolygon* facePtr);

		//// access
		//XPlane* getPlaneBuffer() { return m_pplanes; }
		//const XPlane* getPlaneBuffer() const { return m_pplanes; }

	private:
		MemoryManager() {}
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
}
