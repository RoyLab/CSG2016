#pragma once
#include "macroutil.h"
#include "csgdefs.h"

#include <CGAL\Polyhedron_3.h>
#include <CGAL\Triangle_3.h>
#include <CGAL\HalfedgeDS_vertex_max_base_with_id.h>
#include <CGAL\HalfedgeDS_halfedge_max_base_with_id.h>
#include <CGAL\HalfedgeDS_face_max_base_with_id.h>
#include <list>

// for implementation
#include <CGAL\bounding_box.h>


namespace CSG
{
    struct SharedEdge;
    struct IsectTriangleInfo;

    template <class Refs, class Point>
    struct MyVertex : public CGAL::HalfedgeDS_vertex_max_base_with_id<Refs, Point, unsigned> {
        PointListItrListItr agency = nullptr;
    };

    template <class Refs>
    struct MyHalfedge : public CGAL::HalfedgeDS_halfedge_max_base_with_id<Refs, unsigned> {
        SharedEdge* sharedEdge = nullptr;
        char id012 = -1;
    };

    template <class Refs>
    struct MyFacet : public CGAL::HalfedgeDS_face_max_base_with_id<Refs, CGAL::Tag_false, unsigned> {

        typedef Cube_3 TriBbox;

        CGAL::Triangle_3<K>         triangle;
        Halfedge_handle             edges[3];

        CGAL::Vector_3<K>           normal;
        TriBbox                     box;
        CGAL::Color                 color;

        IsectTriangleInfo*          isectInfo = nullptr;
    };

    struct MyItems {

        template < class Refs, class Traits>
        struct Vertex_wrapper {
            typedef typename Traits::Point_3 Point;
            typedef MyVertex<Refs, Point> Vertex;
        };


        template < class Refs, class Traits>
        struct Halfedge_wrapper {
            typedef MyHalfedge<Refs> Halfedge;
        };

        template <class Refs, class Traits>
        struct Face_wrapper {
            typedef MyFacet<Refs> Face;
        };
    };


    class MyMesh :
        public CGAL::Polyhedron_3<K, MyItems>
    {
        COMMON_PROPERTY(Bbox_3, bbox);
    public:
        void update()
        {
            updateBbox();
            updateFaceNormal();
        }

    private:

        void updateFaceNormal();
        void updateTriangleBbox();

        void updateBbox()
        {
            m_bbox = CGAL::bounding_box(points_begin(), points_end()).bbox();
        }
    };

}