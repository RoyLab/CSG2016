#pragma once
#include "macro_util.h"
#include "csgdefs.h"

#include <CGAL\Polyhedron_3.h>
#include <CGAL\Triangle_3.h>

// for implementation
#include <CGAL\bounding_box.h>


namespace CSG
{

    template <class Refs, class Point>
    struct MyVertex : public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point> {

    };

    template <class Refs>
    struct MyHalfedge : public CGAL::HalfedgeDS_halfedge_base<Refs> {

    };

    template <class Refs>
    struct MyFacet : public CGAL::HalfedgeDS_face_base<Refs> {
        CGAL::Triangle_3<K> triangle;
        //CGAL::Color color;
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