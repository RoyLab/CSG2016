// +-------------------------------------------------------------------------
// | main.cpp
// | 
// | Author: Gilbert Bernstein
// +-------------------------------------------------------------------------
// | COPYRIGHT:
// |    Copyright Gilbert Bernstein 2013
// |    See the included COPYRIGHT file for further details.
// |    
// |    This file is part of the Cork library.
// |
// |    Cork is free software: you can redistribute it and/or modify
// |    it under the terms of the GNU Lesser General Public License as
// |    published by the Free Software Foundation, either version 3 of
// |    the License, or (at your option) any later version.
// |
// |    Cork is distributed in the hope that it will be useful,
// |    but WITHOUT ANY WARRANTY; without even the implied warranty of
// |    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// |    GNU Lesser General Public License for more details.
// |
// |    You should have received a copy 
// |    of the GNU Lesser General Public License
// |    along with Cork.  If not, see <http://www.gnu.org/licenses/>.
// +-------------------------------------------------------------------------

// This file contains a command line program that can be used
// to exercise Cork's functionality without having to write
// any code.

#include "files.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <sstream>
using std::stringstream;
using std::string;

using std::ostream;

#include "cork.h"
#include "corkmesh.h"

#include <xtimer.hpp>
#include <xlogger.h>
#include "stdcsg.h"
#include <vector>
#include <fstream>
using std::vector;

#define XRWY_EXPORTS
#include "cgaleval.h"

namespace
{

    void file2corktrimesh(
        const Files::FileMesh &in, CorkTriMesh *out
    ) {
        out->n_vertices = in.vertices.size();
        out->n_triangles = in.triangles.size();

        out->triangles = new uint[(out->n_triangles) * 3];
        out->vertices = new float[(out->n_vertices) * 3];

        for (uint i = 0; i<out->n_triangles; i++) {
            (out->triangles)[3 * i + 0] = in.triangles[i].a;
            (out->triangles)[3 * i + 1] = in.triangles[i].b;
            (out->triangles)[3 * i + 2] = in.triangles[i].c;
        }

        for (uint i = 0; i<out->n_vertices; i++) {
            (out->vertices)[3 * i + 0] = in.vertices[i].pos.x;
            (out->vertices)[3 * i + 1] = in.vertices[i].pos.y;
            (out->vertices)[3 * i + 2] = in.vertices[i].pos.z;
        }
    }

    void corktrimesh2file(
        CorkTriMesh in, Files::FileMesh &out
    ) {
        out.vertices.resize(in.n_vertices);
        out.triangles.resize(in.n_triangles);

        for (uint i = 0; i<in.n_triangles; i++) {
            out.triangles[i].a = in.triangles[3 * i + 0];
            out.triangles[i].b = in.triangles[3 * i + 1];
            out.triangles[i].c = in.triangles[3 * i + 2];
        }

        for (uint i = 0; i<in.n_vertices; i++) {
            out.vertices[i].pos.x = in.vertices[3 * i + 0];
            out.vertices[i].pos.y = in.vertices[3 * i + 1];
            out.vertices[i].pos.z = in.vertices[3 * i + 2];
        }
    }

    void loadMesh(string filename, CorkTriMesh *out)
    {
        Files::FileMesh filemesh;

        if (Files::readTriMesh(filename, &filemesh) > 0) {
            cerr << "Unable to load in " << filename << endl;
            exit(1);
        }

        file2corktrimesh(filemesh, out);
    }
    void saveMesh(string filename, CorkTriMesh in)
    {
        Files::FileMesh filemesh;

        corktrimesh2file(in, filemesh);

        if (Files::writeTriMesh(filename, &filemesh) > 0) {
            cerr << "Unable to write to " << filename << endl;
            exit(1);
        }
    }

}

namespace csg
{
    class CorkEval : public CsgBinaryEval<CorkMesh>
    {
    public:
        CorkMesh* operator()(CorkMesh* a, CorkMesh* b, CsgOperator op) const
        {
            switch (op)
            {
            case CSG_OP_UNION:
                a->boolUnion(*b);
                break;
            case CSG_OP_INTERSECT:
                a->boolIsct(*b);
                break;
            case CSG_OP_DIFF:
                a->boolDiff(*b);
                break;
            default:
                throw std::exception();
            }
            CorkMesh* result = new CorkMesh(*a);
            return result;
        }
    };
}

using namespace csg;
void corkTriMesh2CorkMesh(
    CorkTriMesh in,
    CorkMesh *mesh_out
);

void corkMesh2CorkTriMesh(
    CorkMesh *mesh_in,
    CorkTriMesh *out
);

extern "C" XRWY_DLL void corkeval(std::vector<std::string>& names, const std::string& expr, 
    const std::string& output, const std::string& log)
{
    XLOG_DEBUG << "\nHere is cgal evaluation.";

    CsgTree<CorkMesh> csg;
    vector<CorkTriMesh*> meshes;
    vector<CorkMesh*> nef_meshes;

    CorkTriMesh* polyheron = nullptr;
    for (string& meshname : names)
    {
        polyheron = new CorkTriMesh;
        loadMesh(meshname,polyheron);
        XLOG_INFO << "Reading off file: " << meshname << " : " << polyheron->n_triangles;
        meshes.push_back(polyheron);
    }

    bool error = false;
    double total = 0;
    int ovn = 0;
    int ofn = 0;

    try
    {
        XTIMER_OBJ.setClock("main");
        XTIMER_OBJ.setClock("step");
        CorkMesh *nefpoly = nullptr;
        for (CorkTriMesh* poly : meshes)
        {
            nefpoly = new CorkMesh;
            corkTriMesh2CorkMesh(*poly, nefpoly);
            nef_meshes.push_back(nefpoly);
        }

        XLOG_INFO << "Convert to  nef polygon: " << XTIMER_OBJ.millisecondsAndReset("step") << " ms";
        csg.createCSGTreeFromExpr(expr, nef_meshes);
        CorkEval eval_obj;
        CorkMesh *nef_result = csg.eval(&eval_obj);

        XLOG_INFO << "eval: " << XTIMER_OBJ.millisecondsAndReset("step") << " ms";
        CorkTriMesh* result = new CorkTriMesh;
        corkMesh2CorkTriMesh(nef_result, result);

        XLOG_INFO << "Convert back: " << XTIMER_OBJ.millisecondsAndReset("step") << " ms";

        total = XTIMER_OBJ.millisecondsAndReset("main");
        XLOG_INFO << "total time: " << total << " ms";

        ofn = result->n_triangles;
        ovn = result->n_vertices;

        delete nef_result;

        for (CorkMesh* nefpoly : nef_meshes)
            delete nefpoly;
        for (CorkTriMesh*poly : meshes)
            delete poly;

        saveMesh(output, *result);
        delete result;
    }
    catch (...)
    {
        error = true;
    }

    if (log.size())
    {
        std::ofstream out(log);
        if (out.is_open())
        {
            out << total << std::endl;
            out << ofn << std::endl;
            out << ovn << std::endl;
            out << error << std::endl;
            out.close();
        }
    }
}








