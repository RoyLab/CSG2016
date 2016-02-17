#pragma once
#include <CGAL\Bbox_3.h>
#include <CGAL\Iso_cuboid_3.h>


namespace myext
{
    CGAL::Bbox_3 enlarge(const CGAL::Bbox_3& bbox, double factor)
    {
        return CGAL::Bbox_3(bbox.xmin() - factor,
            bbox.ymin() - factor,
            bbox.zmin() - factor,
            bbox.xmax() - factor,
            bbox.ymax() - factor,
            bbox.zmax() - factor
            );
    }

    template <class _R>
    bool is_inside_box(const CGAL::Iso_cuboid_3<_R>& boxBig, const CGAL::Iso_cuboid_3<_R>& boxSmall)
    {
        return (boxBig.has_on_bounded_side(boxSmall.min()) &&
            boxBig.has_on_bounded_side(boxSmall.max()));
    }


}