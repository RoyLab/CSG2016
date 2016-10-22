#pragma once
#include <list>

#include <macros.h>

#include "global.h"
#include "preps.h"
#include "intersection.h"

namespace Boolean
{
    class InsctData
    {
    public:
        typedef std::list<PBIRep> Data;

        void refine();
        void isRefined();
        Data& data() { return inscts; }
        const Data& data() const { return inscts; }

    protected:
        bool        refined = false;
        Data        inscts;
    };

    struct Edge: public InsctData
    {

    };

    class IPolygon: public InsctData
    {
    public:

    };

    class Triangle: public IPolygon 
    {

    };

    class Polygon: public IPolygon
    {

    };

	class RegularMesh:
		public ICSGMesh
	{
        typedef cyPointT        PointT;
        typedef IPolygon        FaceT;
        typedef Edge            EdgeT;
	public:
		static RegularMesh* loadFromFile(const char*);
		static RegularMesh* writeFile(const RegularMesh& mesh, const char*);

		bool& inverse();
		const bool& inverse() const;
		size_t id() const;

    protected:
        PointT*     m_vertices;
        FaceT*      m_faces;
        EdgeT*      m_edges;

	};
}
