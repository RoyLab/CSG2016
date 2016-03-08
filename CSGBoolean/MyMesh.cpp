#include <CGAL\bounding_box.h>

#include "MyMesh.h"


namespace CSG
{
    void MyMesh::calcBbox()
    {
        set_bbox(CGAL::bounding_box(points_begin(), points_end()).bbox());
    }

}