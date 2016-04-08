#include <CGAL\bounding_box.h>
#include <cmath>

#include "adaptive.h"
#include "MyMesh.h"


namespace CSG
{
#define NORMALIZE(x, c, s) ((x) * (s) / 2.0 + (c))

    template <class _R>
    inline CGAL::Point_3<_R> normalizePoint(const CGAL::Point_3<_R>& v, double *center, double* scale)
    {
        return CGAL::Point_3<_R>(NORMALIZE(v.x(), center[0], scale[0]),
            NORMALIZE(v.y(), center[1], scale[1]),
            NORMALIZE(v.z(), center[2], scale[2]));
    }

    void MyMesh::calcBbox()
    {
        set_bbox(CGAL::bounding_box(points_begin(), points_end()).bbox());
    }

    void MyMesh::calcTriangles()
    {
        for (auto itr = facets_begin(); itr != facets_end(); itr++)
        {
            auto loop = itr->halfedge();
            K::Point_3 pts[3];
            for (size_t i = 0; i < 3; i++)
            {
                pts[i] = loop->vertex()->point();
                itr->vertices[i] = loop->vertex();
                loop = loop->next();
            }
            itr->data.reset(new UserFData(pts));
            itr->data->sp = Plane_ext<K>(pts[0], pts[1], pts[2]);
        }
    }

    void MyMesh::initEdgeIds()
    {
        for (auto itr = facets_begin(); itr != facets_end(); itr++)
        {
            auto loop = itr->halfedge();
            for (size_t i = 0; i < 3; i++)
            {
                itr->edges[(i+1)%3] = loop;
                loop = loop->next();
            }

            for (size_t i = 0; i < 3; i++)
            {
                assert(itr->edges[i]->vertex() != itr->vertices[i]);
                assert(itr->edges[i]->opposite()->vertex() != itr->vertices[i]);
                loop = loop->next();
            }
        }
    }
    void MyMesh::initIds()
    {
        int count = 0;
        for (auto itr = facets_begin(); itr != facets_end(); itr++)
        {
            itr->id() = count;
            count++;
        }
    }

    void MyMesh::init()
    {
        calcTriangles();
        initEdgeIds();
        initIds();
    }

    /* -1.0 ~ +1.0 */
    void MyMesh::normalize(CGAL::Bbox_3& aabb)
    {
        double center[3], scale[3];
        for (size_t i = 0; i < 3; i++)
        {
            center[i] = (aabb.max(i) + aabb.min(i)) / 2.0;
            scale[i] = aabb.max(i) - aabb.min(i);
        }
        
        for (auto v = points_begin(); v != points_end(); v++)
            *v = normalizePoint(*v, center, scale);

        m_bbox = CGAL::Bbox_3(NORMALIZE(m_bbox.xmin(), center[0], scale[0]),
            NORMALIZE(m_bbox.ymin(), center[1], scale[1]),
            NORMALIZE(m_bbox.zmin(), center[2], scale[2]),
            NORMALIZE(m_bbox.xmax(), center[0], scale[0]),
            NORMALIZE(m_bbox.ymax(), center[1], scale[1]),
            NORMALIZE(m_bbox.zmax(), center[2], scale[2]));
    }

    void MyMesh::denormalize(CGAL::Bbox_3& aabb)
    {
        double center[3], scale[3];
        for (size_t i = 0; i < 3; i++)
        {
            center[i] = (aabb.max(i) + aabb.min(i)) / 2.0;
            scale[i] = aabb.max(i) - aabb.min(i);
        }

        for (auto v = points_begin(); v != points_end(); v++)
        {
            *v = CGAL::Point_3<K>(v->x() * scale[0] / 2.0 + center[0],
                v->y() * scale[1] / 2.0 + center[1],
                v->z() * scale[2] / 2.0 + center[2]);
        }
    }

    void MyMesh::filter(int n)
    {
        const double factor = std::pow(2.0, 18);
        for (auto v = points_begin(); v != points_end(); v++)
        {
            *v = CGAL::Point_3<K>(FP_FILTER(v->x(), factor),
                FP_FILTER(v->y(), factor), FP_FILTER(v->y(), factor));
        }

        m_bbox = CGAL::Bbox_3(FP_FILTER(m_bbox.xmin(), factor),
            FP_FILTER(m_bbox.ymin(), factor),
            FP_FILTER(m_bbox.zmin(), factor),
            FP_FILTER(m_bbox.xmax(), factor),
            FP_FILTER(m_bbox.ymax(), factor),
            FP_FILTER(m_bbox.zmax(), factor));
    }

}