#pragma once
#include <CGAL\Bbox_3.h>
#include <CGAL\Iso_cuboid_3.h>
#include <CGAL\Triangle_3.h>


namespace myext
{
    enum Sign
    {
        UNKOWN = 0,
        INTERSECT_ON_LINE,
        INTERSECT_ON_POINT,
        COPLANAR,
        NOT_INTERSECT
    };

    enum PosTag
    {
        NONE = -1, INNER = 0x00,
        EDGE_0 = 0x01, EDGE_1 = 0x02, EDGE_2 = 0x04,
        VER_0 = 0x08, VER_1 = 0x10, VER_2 = 0x20
    };

    enum BoundedTag
    {
        BT_UNKOWN = 0, BT_NA = 5, BT_ON = 1, BT_INSIDE = 2, BT_OUSIDE = 3
    };

    inline CGAL::Bbox_3 enlarge(const CGAL::Bbox_3& bbox, double factor)
    {
        return CGAL::Bbox_3(bbox.xmin() - factor, bbox.ymin() - factor,
            bbox.zmin() - factor, bbox.xmax() - factor,
            bbox.ymax() - factor, bbox.zmax() - factor);
    }

    template <class _R>
    inline bool is_inside_box(const CGAL::Iso_cuboid_3<_R>& boxBig, const CGAL::Iso_cuboid_3<_R>& boxSmall)
    {
        return (boxBig.has_on_bounded_side(boxSmall.min()) &&
            boxBig.has_on_bounded_side(boxSmall.max()));
    }

    template <class _R>
    inline typename _R::FT dot(const CGAL::Vector_3<_R>& v, const CGAL::Point_3<_R>& p)
    {
        return v.x() * p.x() + v.y() * p.y() + v.z() * p.z();
    }

    template <class _R>
    inline typename _R::FT dot(const CGAL::Point_3<_R>& p, const CGAL::Vector_3<_R>& v)
    {
        return v.x() * p.x() + v.y() * p.y() + v.z() * p.z();
    }

    template <class _R>
    inline int max_coord(const CGAL::Vector_3<_R>& v)
    {
        typedef typename _R::FT FT;
        int index = 0;

        FT max = fabs(LineDir[0]);
        FT b = fabs(LineDir[1]);
        FT c = fabs(LineDir[2]);

        if (b > max) { max = b; index = 1; }
        if (c > max) { max = c; index = 2; }

        return index;
    }
}