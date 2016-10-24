#pragma once
#include <list>
#include <macros.h>

#include "global.h"
#include "preps.h"
#include "intersection.h"
#include "offio.h"

namespace Boolean
{
    struct MyVertex
    {
        int prepId = -1;
        std::list<size_t> edges;
    };

    struct FH
    {
        bool orientation;
        size_t idx;
    };
        
    struct MyEdge
    {
        FH fhs[2];
        FH *extra = 0;
        InsctData* inscts = nullptr;

		~MyEdge() { SAFE_DELETE(inscts); }
    };

    class IPolygon
    {
    public:
		virtual ~IPolygon() {}
        size_t degree() const { return m_degree; }

    protected:
        size_t m_degree = 0;
    };

    class Triangle : public IPolygon
    {
		friend class RegularMesh;
    public:
        bool collide(const CGAL::Iso_cuboid_3<Depick>& cube) const;

    protected:
        Depick::Triangle_3 cgalTri;
        size_t eIds[3];
        size_t vIds[3];

        int spId = -1;
        int bpIds[3] = { -1, -1, -1 };
        InsctData* inscts = nullptr;

		~Triangle() { SAFE_DELETE(inscts); }
	};

    class SubPolygon : public IPolygon
    {

    };

	class RegularMesh:
		public ICSGMesh
	{
        typedef IPolygon        FaceT;
	public:
		static RegularMesh* loadFromFile(const char*);
		static RegularMesh* writeFile(const RegularMesh& mesh, const char*);

		RegularMesh();
		RegularMesh(const OffFile& file); // triangle mesh
		~RegularMesh();

        // csg related
		bool& inverse() { return m_bInverse; }
		const bool& inverse() const { return m_bInverse; }
		size_t id() const { return m_id; }

        // geometry info
        //Bbox_3& bbox();
		void prepareBoolean();

    protected:
        static  MemoryManager* memmgr;

        std::vector<FaceT*>  m_faces;

		int m_id = -1;
		bool m_bInverse = false;
	};
}
