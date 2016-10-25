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
		MyEdge(uint32_t a, uint32_t b)
		{ ends[0] = a; ends[1] = b; }

		uint32_t ends[2];

        FH fhs[2];
		std::list<FH> extrafhs;

        InsctData* inscts = nullptr;

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
		Triangle(size_t i): IPolygon(3, i) {}
        bool collide(const CGAL::Iso_cuboid_3<Depick>& cube) const;
		const CGALTriangle& triangle() const { return cgalTri; }
		XPlane boundingPlane(int i) const { return bPlanes[i]; }

    protected:
		CGALTriangle cgalTri;
        size_t eIds[3];
        size_t vIds[3];

		XPlane sPlane;
		XPlane bPlanes[3];
        InsctData* inscts = nullptr;

		~Triangle() { SAFE_DELETE(inscts); }
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
