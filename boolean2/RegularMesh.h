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
		IPolygon(uint32_t d, uint32_t i):
			m_degree(d), m_id(i){}
		virtual ~IPolygon() {}
		uint32_t degree() const { return m_degree; }
		uint32_t id() const { return m_id; }

    protected:
        const uint32_t m_degree;
		const uint32_t m_id;
    };

    class Triangle : public IPolygon
    {
		friend class RegularMesh;
    public:
		InsctData<FacePBI>* inscts = nullptr;

		Triangle(int i): IPolygon(3, i) {}
        ~Triangle();

		// access
		const CGALTriangle& triangle() const { return cgalTri; }
		XPlane boundingPlane(int i) const { return bPlanes[i]; }
		XPlane supportingPlane() const { return sPlane; }
        cyPointT& point(int i) const;
        MyEdge& edge(int i) const;
        uint32_t edgeId(int i) const { return eIds[i]; }
        uint32_t vertexId(int i) const { return vIds[i]; }

		// search
        uint32_t findVertex(const XPoint& pt, PosTag tag, uint32_t*&);

		// manipulate
		void calcSupportingPlane();
		void calcBoundingPlane();
        template <class Container>
        void addTo(Container& c);

        // state
        bool isAdded4Tess() const { return added; }

    protected:
		CGALTriangle cgalTri;
        uint32_t eIds[3];
        uint32_t vIds[3];

		XPlane sPlane;
		XPlane bPlanes[3];

        bool added = false;
	};

    class SubPolygon : public IPolygon
    {
    public:
        SubPolygon(uint32_t d, uint32_t i = uint32_t(-1)): IPolygon(d, i) {}

        template <class ForwardIterator>
        void constructFromVertexList(const ForwardIterator& a, const ForwardIterator& b);

    protected:
        std::vector<MyEdge::Index> eIds;
        std::vector<MyVertex::Index> vIds;
    };

	class RegularMesh:
		public ICSGMesh
	{
	public:
		typedef IPolygon FaceT;

	protected:
		static  MemoryManager* memmgr;
		std::vector<FaceT*>  m_faces;
		uint32_t m_id;
		bool m_bInverse = false;

	public:
		static RegularMesh* loadFromFile(const char*, uint32_t id);
		static RegularMesh* writeFile(const RegularMesh& mesh, const char*);

		RegularMesh() {}
		RegularMesh(const XR::OffFile& file); // triangle mesh
		~RegularMesh() {}

        // csg related
		bool& inverse() { return m_bInverse; }
		const bool& inverse() const { return m_bInverse; }
		uint32_t id() const { return m_id; }

		// access
		std::vector<FaceT*>& faces() { return m_faces; }
		const std::vector<FaceT*>& faces() const  {return m_faces; }

        // geometry info
        //Bbox_3& bbox();
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

        for (int i = 0; i < degree(); i++)
        {
            v0 = *itr; ++itr;
            v1 = *itr;
            vIds[i] = v0;
            eIds[i] = pMem->getEdgeId(v0, v1, this);
        }
    }
}
