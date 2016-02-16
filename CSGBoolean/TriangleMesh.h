#pragma once
#include <vector>
#include <CGAL\Exact_predicates_inexact_constructions_kernel.h>

namespace CSG
{
    class TriangleMesh
    {
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        typedef CGAL::Point_3<K>                                    Point;
        typedef CGAL::Vector_3<K>                                   Vector;
        typedef float                                               Color[4];
        typedef float                                               Normal[3];

        struct Vertex
        {
            Point pos;

            /** for rendering */
            Normal normal;
        };

        struct Edge
        {
            int fIds[2];
        };

        struct Facet
        {
            int vIds[3];
            int eIds[3];

            /** for rendering */
            Normal normal;
            int colorId;
        };

    protected:
        
        std::vector<Vertex>     m_vertices;
        std::vector<Edge>       m_edges;
        std::vector<Facet>      m_facets;

        int Id;
    };


    bool loadTriangleMeshFromOff(const char* path, TriangleMesh** result, bool checked = true);
}

