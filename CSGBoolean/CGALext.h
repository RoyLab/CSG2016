#pragma once
#include <CGAL\Bbox_3.h>


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
}