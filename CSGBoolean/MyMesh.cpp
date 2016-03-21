#include <CGAL\bounding_box.h>

#include "MyMesh.h"


namespace CSG
{
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
        calcBbox();
        calcTriangles();
        initEdgeIds();
        initIds();
    }

}