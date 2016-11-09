#pragma once
#include <list>
#include <macros.h>

#include "global.h"
#include "preps.h"
#include "offio.h"
#include "xmemory.h"

namespace Boolean
{
    class IPolygon
    {
    public:
        enum TYPE {TRIANGLE, SUBPOLYGON};
        int mark;

    public:
		IPolygon(uint32_t d, uint32_t i, uint32_t meshId):
			m_degree(d), m_id(i), m_meshId(meshId) {}
		virtual ~IPolygon() {}

		uint32_t degree() const { return m_degree; }
		uint32_t id() const { return m_id; }
        uint32_t meshId() const { return m_meshId; }

        virtual bool isValid() const { return true; }
        virtual void getVertices(std::vector<MyVertex::Index>&) const = 0;
        virtual void getEdges(std::vector<MyEdge::Index>&) const;
        virtual TYPE getType() const = 0;

        virtual MyEdge& edge(int i) const;
        virtual uint32_t edgeId(int i) const = 0;
        virtual MyVertex& vertex(int i)const;
        virtual uint32_t vertexId(int i)const = 0;

        XPlane supportingPlane() const { return sPlane; }
        XPlane sPlane;
    protected:
        const uint32_t m_degree;
        const uint32_t m_id;
        const uint32_t m_meshId;
    };

    class Triangle : public IPolygon
    {
		friend class RegularMesh;
    public:
		InsctData<FacePBI>* inscts = nullptr;

		Triangle(uint32_t meshId, uint32_t i): IPolygon(3, i, meshId) {}
        ~Triangle();

		// access
		const CGALTriangle& triangle() const { return cgalTri; }
		XPlane boundingPlane(int i) const { return bPlanes[i]; }
        cyPointT& point(int i) const;
        uint32_t edgeId(int i) const { return eIds[i]; }
        uint32_t vertexId(int i) const { return vIds[i]; }
        void getVertices(std::vector<MyVertex::Index>&) const;
        bool coherentEdge(int whichEdge) const; // 等价于MyEdge的faceOrientation
        MyVertex::Index getTheOtherVertex(MyEdge::Index eId) const;

		// search
        uint32_t findVertex(const XPoint& pt, PosTag tag, uint32_t*&);

		// manipulate
		void calcSupportingPlane();
		void calcBoundingPlane();
        template <class Container>
        void addTo(Container& c);

        // state
        bool isAdded4Tess() const { return added; }
        bool isValid() const { return bIsValid; }
        void invalidate() { bIsValid = false; }
        TYPE getType() const { return TRIANGLE; }

    protected:
		CGALTriangle cgalTri;
        uint32_t eIds[3];
        uint32_t vIds[3];

		XPlane bPlanes[3];

        bool added = false;
        bool bIsValid = true;
	};

    class SubPolygon : public IPolygon
    {
    public:
        SubPolygon(uint32_t meshId, uint32_t d, uint32_t i = uint32_t(-1)):
            IPolygon(d, i, meshId),  eIds(d), vIds(d) {}

        template <class ForwardIterator>
        void constructFromVertexList(const ForwardIterator& a, const ForwardIterator& b);
        void getVertices(std::vector<MyVertex::Index>&) const;
        TYPE getType() const { return SUBPOLYGON; }

        uint32_t edgeId(int i) const { return eIds[i]; }
        uint32_t vertexId(int i) const { return vIds[i]; }

    protected:
        std::vector<MyEdge::Index> eIds;
        std::vector<MyVertex::Index> vIds;
    };

	class RegularMesh:
		public ICSGMesh
	{
	public:
		typedef IPolygon FaceT;
        std::vector<bool> inverseMap;

	protected:
		static  MemoryManager* memmgr;
		std::vector<FaceT*>  m_faces;
		uint32_t m_id;
		bool m_bInverse = false;
        cyPointT m_center, m_scale;

	public:
		static RegularMesh* loadFromFile(const char*, uint32_t id);
		static void writeFile(RegularMesh& mesh, const char*);

		RegularMesh():m_center(0,0,0), m_scale(1,1,1) {}
		RegularMesh(const XR::OffFile& file, uint32_t meshId); // triangle mesh
		~RegularMesh() {}

        // csg related
		bool& inverse() { return m_bInverse; }
		const bool& inverse() const { return m_bInverse; }
		uint32_t id() const { return m_id; }

        /// do not actually change te coordinates
        void invCoords(const cyPointT& c, const cyPointT& s) { m_center = c; m_scale = s; }

		// access
		std::vector<FaceT*>& faces() { return m_faces; }
		const std::vector<FaceT*>& faces() const  {return m_faces; }

        // geometry info
		void prepareBoolean();
	};

    template<class Container>
    inline void Triangle::addTo(Container & c)
    {
        if (isAdded4Tess())
            return;

        c.push_back(this);
        added = true;
    }

    template<class ForwardIterator>
    inline void SubPolygon::constructFromVertexList(const ForwardIterator & a, const ForwardIterator & b)
    {
        ForwardIterator itr = a;
        MyVertex::Index v0, v1;
        auto pMem = MemoryManager::getInstance();

        for (uint32_t i = 0; i < degree(); i++)
        {
            v0 = *itr; ++itr;
            if (itr != b)
                v1 = *itr;
            else
                v1 = *a;

            vIds[i] = v0;
            eIds[i] = pMem->getEdgeId(v0, v1, this);
        }
    }
}
