#include <boost/foreach.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <vector>
#include <fstream>
#include <limits>
#include <string>
#include <fstream>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <boost\log\trivial.hpp>

#define CSG_EXPORTS

#include "bool.h"
#include "MyAlgorithm.h"
#include "MyMesh.h"
#include "adaptive.h"

namespace CSG
{
    extern "C" CSG_API int loadMesh(std::vector<MyMesh*>& meshes, const std::vector<std::string>& names)
    {
        for (int i = 0; i < names.size(); i++)
        {
            std::ifstream file(names[i]);
            assert(file);

            auto newPoly = new MyMesh;
            file >> *newPoly;
            newPoly->init();
            newPoly->Id = i;
            meshes.push_back(newPoly);
        }

        return 0;
    }

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3 Point;
    typedef CGAL::Polyhedron_3<K> Polyhedron;
    typedef MyMesh Polyhedron2;
    using namespace CGAL;
    double max_coordinate(const Polyhedron& poly)
    {
        double max_coord = (std::numeric_limits<double>::min)();
        BOOST_FOREACH(Polyhedron::Vertex_handle v, vertices(poly))
        {
            Point p = v->point();
            max_coord = (std::max)(max_coord, p.x());
            max_coord = (std::max)(max_coord, p.y());
            max_coord = (std::max)(max_coord, p.z());
        }
        return max_coord;
    }

    int main1(const char* filename)
    {
        std::ifstream input(filename);
        Polyhedron poly;
        Polyhedron2 poly2;
        if (!input || !(input >> poly) || poly.empty())
        {
            std::cerr << "Not a valid off file." << std::endl;
            return 1;
        }
        input.seekg(0);
        if (!input || !(input >> poly2) || poly2.empty())
        {
            std::cerr << "Not a valid off file." << std::endl;
            return 1;
        }
        poly2.init();
        input.close();

        CGAL::Side_of_triangle_mesh<Polyhedron, K> inside(poly);
        double size = max_coordinate(poly);
        unsigned int nb_points = 100000;
        std::vector<Point> points;
        points.reserve(nb_points);
        CGAL::Random_points_in_cube_3<Point> gen(size/2);
        for (unsigned int i = 0; i < nb_points; ++i)
            points.push_back(*gen++);
        std::cout << "Test " << nb_points << " random points in cube "
            << "[-" << size << "; " << size << "]" << std::endl;
        int nb_inside = 0;
        int nb_boundary = 0;

        Octree* pOctree = new Octree;
        std::vector<Octree::Node*> intersectLeaves;
        std::vector<MyMesh*> meshList;
        meshList.push_back(&poly2);
        pOctree->build(meshList, &intersectLeaves);

        BOOST_LOG_TRIVIAL(info) << "end of prepration.!\n";

        int code[5] = { 0, 1, -1, 0, 0 };
        for (std::size_t i = 0; i < nb_points; ++i)
        {
            if ((i % 1000) == 0)
            {
                BOOST_LOG_TRIVIAL(info) << "checking " << i << std::endl;
            }
            CGAL::Bounded_side res = inside(points[i]);
            CSG::Relation res2 = CSG::PolyhedralInclusionTest(points[i], pOctree, meshList, 0, false);

            if (res != code[res2])
            {
                std::cout << points[i] << std::endl;
                std::cout << res << '\t' << code[res2] << std::endl;
            }
            if (res == CGAL::ON_BOUNDED_SIDE) { ++nb_inside; }
            if (res == CGAL::ON_BOUNDARY) { ++nb_boundary; }
        }

        delete pOctree;
        std::cerr << "Total query size: " << points.size() << std::endl;
        std::cerr << "  " << nb_inside << " points inside " << std::endl;
        std::cerr << "  " << nb_boundary << " points on boundary " << std::endl;
        std::cerr << "  " << points.size() - nb_inside - nb_boundary << " points outside " << std::endl;
        return 0;
    }

    extern "C" CSG_API void test()
    {
        //main1("../../models/horse.off");
        GS::exactinit();
        std::vector<std::string> names;
        names.push_back("../../models/ball1.off");
        names.push_back("../../models/ball2.off");

        std::string expr("0+1");

        std::vector<MyMesh*> meshList;
        loadMesh(meshList, names);

        MyAlgorithm *alg = new MyAlgorithm;
        alg->solve(expr, meshList);

        auto res = alg->getResultMesh();
        std::ofstream output("../../models/ball-res.off");
        output << *res.get();
        output.close();

        for (auto mesh : meshList)
        {
            //std::cout << mesh->size_of_facets() << std::endl;
            SAFE_DELETE(mesh);
        }

        SAFE_DELETE(alg);
    }
}