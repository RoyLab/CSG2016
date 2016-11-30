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

using std::string;
using std::vector;

namespace boolean
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

using namespace boolean;

extern "C" XRWY_DLL void cgaleval(std::vector<std::string>& names, std::string& expr, const std::string& output)
{
    XLOG_DEBUG << "Here is cgal evaluation.";

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
        meshes.push_back(polyheron);
    }

    Nef_polyhedron *nefpoly = nullptr;
    for (Polyhedron* poly : meshes)
    {
        nefpoly = new Nef_polyhedron(*poly);
        nef_meshes.push_back(nefpoly);
    }

    csg.createCSGTreeFromExpr(expr, nef_meshes);
    CgalEval eval_obj;
    Nef_polyhedron *nef_result = csg.eval(&eval_obj);

    Polyhedron* result = new Polyhedron;
    nef_result->convert_to_polyhedron(*result);
    delete nef_result;

    std::ofstream output_file(output);
    output_file << *result;
    output_file.close();

    delete result;
}

extern "C" XRWY_DLL void cgaleval_test()
{
}
