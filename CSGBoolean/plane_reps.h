#pragma once
#include "CGALext.h"
#include "macroutil.h"
#include <CGAL\intersections.h>

namespace CSG
{
    template <class _R>
    class PBPoint
    {
        DECLARE_CGAL_KERNEL_CLASS
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

    private:
        Plane planes[3];
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
        COMMON_PROPERTY(Direction, normal);
        COMMON_PROPERTY(Plane, sp);
        Plane m_bps[3];
    };
}