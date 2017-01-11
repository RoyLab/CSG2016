#include "precompile.h"
#define XRWY_EXPORTS
#include "boolean.h"
#include "preps.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/bounding_box.h>
#include <boost\log\trivial.hpp>
#include <vector>
#include <fstream>
#include <limits>
#include <boost/foreach.hpp>

#include "RegularMesh.h"
#include "Octree.h"

using namespace Boolean;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron;
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

std::string print_text(int i)
{
    switch (i)
    {
    case -1:
        return "out";
    case 0:
        return "on";
    case 1:
        return "in";
    default:
        throw 1;
    }
}

XRWY_DLL void test1()
{
    const char* filename = "D:\\Codes\\Boolean2016\\models\\bunny.off";
    //const char* filename = "D:\\Codes\\Boolean2016\\models\\ball1.off";

    Polyhedron poly;
    std::ifstream input(filename);
    if (!input || !(input >> poly) || poly.empty())
    {
        std::cerr << "Not a valid off file." << std::endl;
        return;
    }
    input.close();

    RegularMesh* poly2 = RegularMesh::loadFromFile(filename, 0);
    if (poly2->faces().empty())
    {
        std::cerr << "Not a valid off file." << std::endl;
        return;
    }

    CGAL::Side_of_triangle_mesh<Polyhedron, K> inside(poly);
    double size = max_coordinate(poly);
    unsigned int nb_points = 100000;

    // sample mode
    //nb_points = 1;
    //K::Point_3 sample(-3.02041, 3.50603, 3.58144);

    std::vector<Point> points;
    points.reserve(nb_points);
    CGAL::Random_points_in_cube_3<Point> gen(size / 2);
    for (unsigned int i = 0; i < nb_points; ++i)
        points.push_back(*gen++);
    std::cout << "Test " << nb_points << " random points in cube "
        << "[-" << size << "; " << size << "]" << std::endl;
    int nb_inside = 0;
    int nb_boundary = 0;

    Octree* pOctree = new Octree;
    std::vector<Octree::Node*> intersectLeaves;
    std::vector<RegularMesh*> meshList;
    meshList.push_back(poly2);
    CGAL::Bbox_3 bbox = CGAL::bounding_box(poly.points_begin(), poly.points_end()).bbox();

    pOctree->build(meshList,
        bbox, 
        true,
        &intersectLeaves);

    BOOST_LOG_TRIVIAL(info) << "end of prepration.!";

    int code[5] = { 0, 1, -1, 0, 0 };
    for (std::size_t i = 0; i < nb_points; ++i)
    {
        if ((i % 1000) == 0)
        {
            BOOST_LOG_TRIVIAL(info) << "checking " << i;
        }

        auto thiz = points[i];
        //auto thiz = sample;
        int cross;
        CGAL::Bounded_side res = inside(thiz);
        Relation res2 = PolyhedralInclusionTest(
            cyPointT(thiz.x(), thiz.y(), thiz.z()),
            pOctree, meshList, 0, false, &cross);

        if (res != code[res2])
        {
            std::cout << thiz << std::endl;
            std::cout << "cgal: " << print_text(res) 
                << '\t' << "my: " << print_text(code[res2]) 
                << " cross: "<< cross << std::endl;
        }
        if (res == CGAL::ON_BOUNDED_SIDE) { ++nb_inside; }
        if (res == CGAL::ON_BOUNDARY) { ++nb_boundary; }
    }

    delete pOctree;
    std::cerr << "Total query size: " << points.size() << std::endl;
    std::cerr << "  " << nb_inside << " points inside " << std::endl;
    std::cerr << "  " << nb_boundary << " points on boundary " << std::endl;
    std::cerr << "  " << points.size() - nb_inside - nb_boundary << " points outside " << std::endl;


    //std::ifstream input(filename);
    //Polyhedron poly;
    //if (!input || !(input >> poly) || poly.empty())
    //{
    //    std::cerr << "Not a valid off file." << std::endl;
    //    return 1;
    //}
    //CGAL::Side_of_triangle_mesh<Polyhedron, K> inside(poly);
    //double size = max_coordinate(poly);
    //unsigned int nb_points = 100;
    //std::vector<Point> points;
    //points.reserve(nb_points);
    //CGAL::Random_points_in_cube_3<Point> gen(size);
    //for (unsigned int i = 0; i < nb_points; ++i)
    //    points.push_back(*gen++);
    //std::cout << "Test " << nb_points << " random points in cube "
    //    << "[-" << size << "; " << size << "]" << std::endl;
    //int nb_inside = 0;
    //int nb_boundary = 0;
    //for (std::size_t i = 0; i < nb_points; ++i)
    //{
    //    CGAL::Bounded_side res = inside(points[i]);
    //    if (res == CGAL::ON_BOUNDED_SIDE) { ++nb_inside; }
    //    if (res == CGAL::ON_BOUNDARY) { ++nb_boundary; }
    //}
    //std::cerr << "Total query size: " << points.size() << std::endl;
    //std::cerr << "  " << nb_inside << " points inside " << std::endl;
    //std::cerr << "  " << nb_boundary << " points on boundary " << std::endl;
    //std::cerr << "  " << points.size() - nb_inside - nb_boundary << " points outside " << std::endl;
    return;
}


//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef K::Point_3 Point;
//typedef CGAL::Polyhedron_3<K> Polyhedron;
//typedef MyMesh Polyhedron2;
//using namespace CGAL;
//double max_coordinate(const Polyhedron& poly)
//{
//    double max_coord = (std::numeric_limits<double>::min)();
//    BOOST_FOREACH(Polyhedron::Vertex_handle v, vertices(poly))
//    {
//        Point p = v->point();
//        max_coord = (std::max)(max_coord, p.x());
//        max_coord = (std::max)(max_coord, p.y());
//        max_coord = (std::max)(max_coord, p.z());
//    }
//    return max_coord;
//}

//int main1(const char* filename)
//{
//    std::ifstream input(filename);
//    Polyhedron poly;
//    Polyhedron2 poly2;
//    if (!input || !(input >> poly) || poly.empty())
//    {
//        std::cerr << "Not a valid off file." << std::endl;
//        return 1;
//    }
//    input.seekg(0);
//    if (!input || !(input >> poly2) || poly2.empty())
//    {
//        std::cerr << "Not a valid off file." << std::endl;
//        return 1;
//    }
//    poly2.init();
//    input.close();

//    CGAL::Side_of_triangle_mesh<Polyhedron, K> inside(poly);
//    double size = max_coordinate(poly);
//    unsigned int nb_points = 100000;
//    std::vector<Point> points;
//    points.reserve(nb_points);
//    CGAL::Random_points_in_cube_3<Point> gen(size/2);
//    for (unsigned int i = 0; i < nb_points; ++i)
//        points.push_back(*gen++);
//    std::cout << "Test " << nb_points << " random points in cube "
//        << "[-" << size << "; " << size << "]" << std::endl;
//    int nb_inside = 0;
//    int nb_boundary = 0;

//    Octree* pOctree = new Octree;
//    std::vector<Octree::Node*> intersectLeaves;
//    std::vector<MyMesh*> meshList;
//    meshList.push_back(&poly2);
//    pOctree->build(meshList, &intersectLeaves);

//    BOOST_LOG_TRIVIAL(info) << "end of prepration.!\n";

//    int code[5] = { 0, 1, -1, 0, 0 };
//    for (std::size_t i = 0; i < nb_points; ++i)
//    {
//        if ((i % 1000) == 0)
//        {
//            BOOST_LOG_TRIVIAL(info) << "checking " << i << std::endl;
//        }
//        CGAL::Bounded_side res = inside(points[i]);
//        CSG::Relation res2 = CSG::PolyhedralInclusionTest(points[i], pOctree, meshList, 0, false);

//        if (res != code[res2])
//        {
//            std::cout << points[i] << std::endl;
//            std::cout << res << '\t' << code[res2] << std::endl;
//        }
//        if (res == CGAL::ON_BOUNDED_SIDE) { ++nb_inside; }
//        if (res == CGAL::ON_BOUNDARY) { ++nb_boundary; }
//    }

//    delete pOctree;
//    std::cerr << "Total query size: " << points.size() << std::endl;
//    std::cerr << "  " << nb_inside << " points inside " << std::endl;
//    std::cerr << "  " << nb_boundary << " points on boundary " << std::endl;
//    std::cerr << "  " << points.size() - nb_inside - nb_boundary << " points outside " << std::endl;
//    return 0;
//}

