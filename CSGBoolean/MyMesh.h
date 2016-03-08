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
    struct VertexDelegate;
    struct SharedEdge;
    struct IsectTriangleInfo;

    struct ContextNode
    {
        unsigned    meshId;
        std::vector<MyMesh::Face_handle> facets;
    };

    typedef std::vector<ContextNode*>  VertexContext;

    enum Mark { UNVISITED, SEEDED, VISITED };

    template <class Refs, class Point>
    struct MyVertex : public CGAL::HalfedgeDS_vertex_max_base_with_id<Refs, Point, unsigned> {
        union 
        {
            VertexDelegate* shared = nullptr;
            VertexContext* ctx = nullptr;
        };
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
        CGAL::Vector_3<K>           normal;

        CGAL::Color                           color;
        
        Halfedge_handle                 edges[3];
        TriBbox                                     box;

        IsectTriangleInfo*              isectInfo = nullptr;

        // for flood filling
        Mark            mark = UNVISITED;
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

    typedef MyMesh::Vertex_iterator     PointListItr;
    typedef std::list<PointListItr>     PointListItrList;
    typedef PointListItrList::iterator  PointListItrListItr;

    struct VertexDelegate
    {
        PointListItrListItr    agency;
        PointListItr           location = nullptr;
    };
    
    struct SharedEdge
    {
        std::vector<PointListItrListItr>  innerPoints;
    };

    struct IsectTriangleInfo
    {
        std::vector< >  innerPoints;
        std::vector<MyMesh::Face_handle> coplanars;
    };
}