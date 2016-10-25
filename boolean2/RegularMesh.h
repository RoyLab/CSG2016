#pragma once
#include <list>
#include <macros.h>

#include "global.h"
#include "preps.h"
#include "intersection.h"
#include "offio.h"

namespace Boolean
{
	typedef Depick CGALKernel;
	typedef typename CGALKernel::Triangle_3 CGALTriangle;
	typedef typename CGALKernel::Point_3 CGALPoint;

    struct MyVertex
    {
		int rep; // positive is v-base, negative is p-base
        std::list<uint32_t> edges;

		bool findEdge(uint32_t other, uint32_t* result = nullptr) const;
		Oriented_side orientation(const XPlane&) const;
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

		MyEdge(uint32_t a, uint32_t b)
		{ ends[0] = a; ends[1] = b; }

		uint32_t ends[2];

        FH fhs[2];
		std::list<FH> extrafhs;
		InsctData<EdgePBI>* inscts = nullptr;

		~MyEdge() { SAFE_DELETE(inscts); }
		void addAjacentFace(uint32_t s, uint32_t e, IPolygon* fPtr);
    };

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
		~Triangle() { SAFE_DELETE(inscts); }

		// access
		const CGALTriangle& triangle() const { return cgalTri; }
		XPlane boundingPlane(int i) const { return bPlanes[i]; }
		XPlane supportingPlane() const { return sPlane; }
		CGALPoint& point(int i) const;
		MyEdge& edge(int i) const;
		uint32_t edgeId(int i) const;

		// search
		int findVertex(const XPoint& pt, PosTag tag) const;

		// manipulate
		void calcSupportingPlane();
		void calcBoundingPlane();

    protected:
		CGALTriangle cgalTri;
        size_t eIds[3];
        size_t vIds[3];

		XPlane sPlane;
		XPlane bPlanes[3];
	};

    class SubPolygon : public IPolygon
    {

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
}
