#pragma once
#include <CGAL\Polyhedron_3.h>
#include <CGAL\HalfedgeDS_vertex_max_base_with_id.h>
#include <CGAL\HalfedgeDS_halfedge_max_base_with_id.h>
#include <CGAL\HalfedgeDS_face_max_base_with_id.h>
#include <set>
#include <map>

#include <boost\shared_ptr.hpp>

#include "macroutil.h"
#include "csgdefs.h"
#include "plane_reps.h"

namespace CSG
{
    struct ItstTriangle;
    struct VEntity;
    class ItstAlg;

    typedef std::list<VEntity*>     VEntities;
    typedef VEntities::iterator     VProxy;
    typedef std::list<VProxy>       VProxies;

    class VProxyItr :
        public std::list<VProxy>::iterator
    {
    public:
        VProxyItr(){}
        VProxyItr(std::list<VProxy>::iterator itr):
            std::list<VProxy>::iterator(itr){}

        VEntity* pointer()
        {
            return ***this;
        }

        const VEntity* pointer() const
        {
            return ***this;
        }
    };


    struct UserVData
    {
        UserVData(VProxyItr p) :proxy(new VProxyItr(p)){}

        VProxyItr* proxy = nullptr;
        ~UserVData(){ SAFE_DELETE(proxy);}
    };

    struct UserEData
    {
        std::vector<VProxyItr>  vertices;
    };

    struct UserFData
    {
        typedef Cube_3 TriBbox;

        UserFData(K::Point_3 pts[]) :
            triangle(pts[0], pts[1], pts[2])
        { 
            bbox = triangle.bbox();
            sp = Plane_ext<K>(pts[0], pts[1], pts[2]);
        }

        CGAL::Triangle_3<K>     triangle;
        Plane_ext<K>            sp;
        PBTriangle<K>*          planeRep = nullptr;
        TriBbox                 bbox;

        ItstTriangle*           itstTri = nullptr;
    };

    template <class Refs, class Point>
    struct MyVertex : public CGAL::HalfedgeDS_vertex_max_base_with_id<Refs, Point, unsigned> {
        MyVertex(){}
        MyVertex(const Point& p):CGAL::HalfedgeDS_vertex_max_base_with_id<Refs, Point, unsigned>(p){}
        boost::shared_ptr<UserVData> data;
        int idx = -1;
    };

    template <class Refs>
    struct MyHalfedge : public CGAL::HalfedgeDS_halfedge_max_base_with_id<Refs, unsigned> {
        boost::shared_ptr<UserEData> data;
        int id = -1;
    };

    template <class Refs>
    struct MyFacet : public CGAL::HalfedgeDS_face_max_base_with_id<Refs, CGAL::Tag_false, unsigned> {
        boost::shared_ptr<UserFData>    data;

        Halfedge_handle                 edges[3];
        Vertex_handle                   vertices[3];
        CGAL::Color                     color;
        int                             mark = MARK_BEGIN;

        bool isSimple() const;
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
        void init();
        void calcBbox();
        void normalize(CGAL::Bbox_3& aabb);
        void denormalize(CGAL::Bbox_3& aabb);
        void filter(int n = 18);

        bool bInverse = false;
        int Id = -1;

    private:
        void initDebugInf();
        void calcTriangles();
        void initEdgeIds();
        void initIds();
    };

#define DEFINE_HANDLES\
    typedef MyMesh::Vertex_handle VH;\
    typedef MyMesh::Halfedge_handle EH;\
    typedef MyMesh::Face_handle FH

    DEFINE_HANDLES;

    /* 注意确保pts顺序顺着方向n0 x plane.sp */
    struct ItstLine
    {
        // tag 和 idx 共同确定了独一无二的面内点坐标
        struct
        {
            //PosTag tag = NONE;
            //int idx = -1;
            VProxyItr vertex;
        } pts[2];

        FH plane;
        int id = -1;

        bool check(const Plane_ext<K>& sp);
    };

    typedef std::list<ItstLine> ItstLineList;
    typedef std::map<int, ItstLineList> ItstLineLists;

    struct ItstTriangle
    {
        ItstLineLists           isectLines;
        ItstLineList            unifiedLines;
        std::vector<VProxyItr>  inVertices;
        std::set<int>           meshIds;


        ItstTriangle(MyMesh::Face_handle fh){}
    };

    template <class Refs>
    bool MyFacet<Refs>::isSimple() const
    {
        return !data->itstTri || data->itstTri->meshIds.empty();
    }

    inline bool getSharedEdge(FH f0, FH f1, EH& eh)
    {
        for (size_t i = 0; i < 3; i++)
        {
            if (f0->edges[i]->opposite()->facet() == f1)
            {
                eh = f0->edges[i];
                return true;
            }
        }
        return false;
    }

#define DECLARE_MYMESH_TYPES\
    typedef MyMesh::Face_handle FH;\
    typedef MyMesh::Halfedge_handle EH;\
    typedef MyMesh::Vertex_handle VH
}