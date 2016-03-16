#pragma once
#include "macroutil.h"
#include "csgdefs.h"

#include <CGAL\Polyhedron_3.h>
#include <CGAL\HalfedgeDS_vertex_max_base_with_id.h>
#include <CGAL\HalfedgeDS_halfedge_max_base_with_id.h>
#include <CGAL\HalfedgeDS_face_max_base_with_id.h>

#include <boost\smart_ptr.hpp>


namespace CSG
{
    struct UserVData;
    struct UserEData;
    struct UserFData;

#define MARK_BEGIN 0xff

    template <class Refs, class Point>
    struct MyVertex : public CGAL::HalfedgeDS_vertex_max_base_with_id<Refs, Point, unsigned> {
        boost::scoped_ptr<UserVData> data;
    };

    template <class Refs>
    struct MyHalfedge : public CGAL::HalfedgeDS_halfedge_max_base_with_id<Refs, unsigned> {
        boost::scoped_ptr<UserEData> data;
        int id = -1;
    };

    template <class Refs>
    struct MyFacet : public CGAL::HalfedgeDS_face_max_base_with_id<Refs, CGAL::Tag_false, unsigned> {
        boost::scoped_ptr<UserFData>    data;

        Halfedge_handle                 edges[3];
        CGAL::Vector_3<K>               normal;
        CGAL::Color                     color;
        int                             mark = -1;
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
        Cube_3 get_bbox_cube() const { return Cube_3(m_bbox); }
        void calcBbox();
    };
}