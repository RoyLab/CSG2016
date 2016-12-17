#define CGAL_NO_ASSERTIONS
#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Cartesian<CGAL::Exact_rational>  EiKernel;
typedef CGAL::Polyhedron_3<EiKernel> Polyhedron;
typedef CGAL::Nef_polyhedron_3<EiKernel>  Nef_polyhedron;


#include <fstream>
#include <string>
#include <xlogger.h>
#include "stdcsg.h"

#define XRWY_EXPORTS
#include "cgaleval.h"

#define XTIMER_INSTANCE
#include <xtimer.hpp>

using std::string;
using std::vector;

namespace csg
{
    class CgalEval : public CsgBinaryEval<Nef_polyhedron>
    {
    public:
        Nef_polyhedron* operator()(Nef_polyhedron* a, Nef_polyhedron* b, CsgOperator op) const
        {
            Nef_polyhedron* result = new Nef_polyhedron;
            switch (op)
            {
            case CSG_OP_UNION:
                *result = *a + *b;
                //std::cout << "Number of halfedges: " << result->number_of_halfedges() << std::endl;
                break;
            case CSG_OP_INTERSECT:
                *result = *a * *b;
                break;
            case CSG_OP_DIFF:
                *result = *a - *b;
            default:
                break;
            }
            return result;
        }
    };
}

using namespace csg;

extern "C" XRWY_DLL void cgaleval(std::vector<std::string>& names, std::string& expr, const std::string& output)
{
    XLOG_DEBUG << "\nHere is cgal evaluation.";

    CsgTree<Nef_polyhedron> csg;
    vector<Polyhedron*> meshes;
    vector<Nef_polyhedron*> nef_meshes;

    Polyhedron* polyheron = nullptr;
    for (string& meshname : names)
    {
        std::ifstream file(meshname);
        polyheron = new Polyhedron;
        file >> *polyheron;
        file.close();
        XLOG_INFO << "Reading off file: " << meshname << " : " << polyheron->size_of_facets();
        meshes.push_back(polyheron);
    }

    XTIMER_OBJ.setClock("main");
    XTIMER_OBJ.setClock("step");
    Nef_polyhedron *nefpoly = nullptr;
    for (Polyhedron* poly : meshes)
    {
        nefpoly = new Nef_polyhedron(*poly);
        nef_meshes.push_back(nefpoly);
    }

    XLOG_INFO << "Convert to  nef polygon: " << XTIMER_OBJ.millisecondsAndReset("step") << " ms";
    csg.createCSGTreeFromExpr(expr, nef_meshes);
    CgalEval eval_obj;
    Nef_polyhedron *nef_result = csg.eval(&eval_obj);

    XLOG_INFO << "eval: " << XTIMER_OBJ.millisecondsAndReset("step") << " ms";
    Polyhedron* result = new Polyhedron;
    nef_result->convert_to_polyhedron(*result);
    delete nef_result;

    XLOG_INFO << "Convert back: " << XTIMER_OBJ.millisecondsAndReset("step") << " ms";
    XLOG_INFO << "total time: " << XTIMER_OBJ.millisecondsAndReset("main") << " ms";

    for (Nef_polyhedron* nefpoly : nef_meshes)
        delete nefpoly;
    for (Polyhedron*poly : meshes)
        delete poly;

    std::ofstream output_file(output);
    output_file << *result;
    output_file.close();

    delete result;
}

extern "C" XRWY_DLL void cgaleval_test()
{
}
