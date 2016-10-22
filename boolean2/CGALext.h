#pragma once
#include <CGAL\Bbox_3.h>
#include <CGAL\Iso_cuboid_3.h>
#include <CGAL\Triangle_3.h>


namespace Boolean
{
    enum BoundedTag
    {
        BT_UNKOWN = 0, BT_NA = 5, BT_ON = 1, BT_INSIDE = 2, BT_OUSIDE = 3
    };

    inline CGAL::Bbox_3 enlarge(const CGAL::Bbox_3& bbox, double space)
    {
        return CGAL::Bbox_3(bbox.xmin() - space, bbox.ymin() - space,
            bbox.zmin() - space, bbox.xmax() + space,
            bbox.ymax() + space, bbox.zmax() + space);
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

#define DECLARE_CGAL_KERNEL_CLASS\
    typedef typename _R::Plane_3 Plane;\
    typedef typename _R::Point_3 Point;\
    typedef typename _R::Vector_3 Vector;\
    typedef typename _R::Triangle_3 Triangle;\
    typedef typename _R::Direction_3 Direction;

}