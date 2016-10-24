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
        int prepId = -1;
        std::list<size_t> edges;
    };

    struct FH
    {
        bool orientation;
        IPolygon* ptr = 0;
    };
        
    struct MyEdge
    {
		size_t points[2];
        FH fhs[2];
		std::list<FH> extrafhs;
        InsctData* inscts = nullptr;

		~MyEdge() { SAFE_DELETE(inscts); }
		size_t addAjacentFace(int s, int e);
    };

    class IPolygon
    {
    public:
		IPolygon(int d, size_t i):
			m_degree(d), m_id(i){}
		virtual ~IPolygon() {}
		size_t degree() const { return m_degree; }
		size_t id() const { return m_id; }

    protected:
        const size_t m_degree;
		const size_t m_id;
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
		bool m_bInverse;

	public:
		static RegularMesh* loadFromFile(const char*);
		static RegularMesh* writeFile(const RegularMesh& mesh, const char*);

		RegularMesh();
		RegularMesh(const OffFile& file); // triangle mesh
		~RegularMesh();

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
