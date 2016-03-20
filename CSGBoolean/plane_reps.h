#pragma once
#include "CGALext.h"
#include "macroutil.h"
#include "adaptive.h"
#include <CGAL\intersections.h>

namespace CSG
{
    template <class _R>
    class Plane_ext :
        CGAL::Plane_3 <_R>
    {
    public:
        Plane_ext(const CGAL::Plane_3& p)
        {
            memset(this, &p, sizeof(Plane_ext));
        }

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
        }
    };

    template <class _R>
    class PBPoint
    {
        DECLARE_CGAL_KERNEL_CLASS
        typedef Plane_ext<_R> PlaneExt;
    public:
        PBPoint(const Plane& p, const Plane& q, const Plane& r)
        {
            bps[0] = p;
            bps[1] = q;
            bps[2] = r;
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
                    mat[i][j] = bps[i][j];

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 4; j++)
                    mat[3][j] = p.bps[i][j];

                if (GS::adaptiveDet4x4Sign(mat) != 0)
                    return false;
            }
            return true;
        }

    private:
        PlaneExt planes[3];
    };

    template <class _R>
    class PBTriangle
    {
        DECLARE_CGAL_KERNEL_CLASS
    public:
        PBTriangle(const Point& p, const Point& q, const Point& r)
        {
            Vector e1 = q - p;
            Vector e2 = r - p;

            m_normal = CGAL::cross_product(e1, e2);
            sp = Plane(p, normal);

            // TO-DO
        }



    private:
        COMMON_PROPERTY(Plane, sp);
        Plane m_bps[3];
    };
}