#pragma once
#include "CGALext.h"
#include "macroutil.h"
#include "adaptive.h"
#include <CGAL\intersections.h>

namespace CSG
{
    template <class _R>
    class Plane_ext :
        public CGAL::Plane_3<_R>
    {
        DECLARE_CGAL_KERNEL_CLASS
    public:
        Plane_ext(){}
        Plane_ext(const CGAL::Plane_3<_R>& p)
        { memcpy(this, &p, sizeof(Plane_ext)); }

        Plane_ext(const Point& p, const Point& q, const Point& r):
            CGAL::Plane_3<_R>(p, CGAL::cross_product(q - p, r - p))
        {}

        typename _R::FT operator[](int i) const
        {
            assert(i < 4 && i >= 0);
            switch (i)
            {
            case 0: return a();
            case 1: return b();
            case 2: return c();
            case 3: return d();
            }
            ReportError();
            return 0;
        }
    };

    template <class _R>
    class PBPoint
    {
        DECLARE_CGAL_KERNEL_CLASS
        typedef Plane_ext<_R> PlaneExt;
    public:
        PBPoint(){}
        PBPoint(const Plane& p, const Plane& q, const Plane& r, bool keepPositive = true)
        {
            planes[0] = p;
            planes[1] = q;
            planes[2] = r;

            if (keepPositive)
            {
                if (CGAL::determinant(p.orthogonal_vector(),
                    q.orthogonal_vector(), r.orthogonal_vector()) < 0.f)
                    planes[2] = r.opposite();
            }

            assert(CGAL::determinant(p.orthogonal_vector(),
                q.orthogonal_vector(), r.orthogonal_vector()) > 0.f);
        }

        Point estimate_coords() const
        {
            return CGAL::intersection(planes[0], planes[1], planes[2]);
        }

        bool operator==(const PBPoint& p) const
        {
            double4x4 mat;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 4; j++)
                    mat[i][j] = p.planes[i][j];

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 4; j++)
                    mat[3][j] = p.planes[i][j];

                if (GS::adaptiveDet4x4Sign(mat) != 0)
                    return false;
            }
            return true;
        }

    private:
        PlaneExt planes[3];
#ifdef _DEBUG
    public:
        void computeCoord()
        {
            auto result = CGAL::intersection(planes[0], planes[1], planes[2]);
            const Point* p = boost::get<Point>(&*result);
            memcpy(&coord, p, sizeof(double3));
        }

        double3 coord;
#endif
    };

    //typedef K _R;
    template <class _R>
    class PBTriangle
    {
        DECLARE_CGAL_KERNEL_CLASS
        typedef Plane_ext<_R> PlaneExt;
    public:
        PBTriangle(){}

        PBTriangle(const Triangle& t)
        {
            new(this)PBTriangle(t.vertex(0), t.vertex(1), t.vertex(2));
        }

        PBTriangle(const Point& p, const Point& q, const Point& r)
        {
            m_sp = PlaneExt(p, q, r);

            Vector normal = m_sp.orthogonal_vector();
            normal = normal / sqrt(normal.squared_length());

            m_bps[0] = PlaneExt(q, r, q + normal);
            m_bps[1] = PlaneExt(r, p, r + normal);
            m_bps[2] = PlaneExt(p, q, p + normal);
        }

    private:
        COMMON_PROPERTY(PlaneExt, sp);
        PlaneExt m_bps[3];
    };

    template <class Plane>
    double orientation(const Plane& p, const Plane& q, const Plane& r, const Plane& s)
    {
        return CGAL::determinant(p.a(), p.b(), p.c(), p.d(),
            q.a(), q.b(), q.c(), q.d(),
            r.a(), r.b(), r.c(), r.d(),
            s.a(), s.b(), s.c(), s.d());
    }
}