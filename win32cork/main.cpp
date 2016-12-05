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
#define XTIMER_INSTANCE
#include <xtimer.hpp>
#include <xlogger.h>
#include <cmdline.h>
#include "stdcsg.h"
#include <vector>
#include <fstream>
using std::vector;

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

void corkeval(std::vector<std::string>& names, std::string& expr, const std::string& output)
{
    XLOG_DEBUG << "\nHere is cork evaluation.";

    CsgTree<CorkMesh> csg;
    vector<CorkTriMesh*> meshes;
    vector<CorkMesh*> nef_meshes;

    CorkTriMesh* polyheron = nullptr;
    for (string& meshname : names)
    {
        polyheron = new CorkTriMesh;
        loadMesh(meshname, polyheron);
        XLOG_INFO << "Reading off file: " << meshname << " : " << polyheron->n_triangles;
        meshes.push_back(polyheron);
    }

    XTIMER_OBJ.setClock("main");
    XTIMER_OBJ.setClock("step");
    CorkMesh *nefpoly = nullptr;
    for (CorkTriMesh* poly : meshes)
    {
        nefpoly = new CorkMesh;
        corkTriMesh2CorkMesh(*poly, nefpoly);
        nef_meshes.push_back(nefpoly);
    }

    XLOG_INFO << "Convert to  cork polygon: " << XTIMER_OBJ.millisecondsAndReset("step") << " ms";
    csg.createCSGTreeFromExpr(expr, nef_meshes);
    CorkEval eval_obj;
    CorkMesh *nef_result = csg.eval(&eval_obj);

    XLOG_INFO << "eval: " << XTIMER_OBJ.millisecondsAndReset("step") << " ms";
    CorkTriMesh* result = new CorkTriMesh;
    corkMesh2CorkTriMesh(nef_result, result);
    delete nef_result;

    XLOG_INFO << "Convert back: " << XTIMER_OBJ.millisecondsAndReset("step") << " ms";
    XLOG_INFO << "total time: " << XTIMER_OBJ.millisecondsAndReset("main") << " ms";

    for (CorkMesh* nefpoly : nef_meshes)
        delete nefpoly;
    for (CorkTriMesh*poly : meshes)
        delete poly;

    saveMesh(output, *result);

    delete result;
}



std::stringstream getValidStream(std::ifstream& file)
{
    std::string line;
    while (file)
    {
        std::getline(file, line);
        if (line.size() && line[0] != ';')
            return (std::stringstream(line));
    }
    throw std::exception();
}

int main(int argc, char* argv[])
{
    cmdline::parser cmd_parser;
    cmd_parser.add<std::string>("script", 's', "the script to run", false, "./mycsg.ini");

    cmd_parser.parse_check(argc, argv);
    std::string config_filename = cmd_parser.get<std::string>("script");

    std::string expr;
    std::vector<std::string> names;

    if (false && argc == 1)
    {
        expr = "0+1";
        names.clear();
        names.push_back("D:/Codes/Boolean2016/exps/data/cmp_meshworks/buddha.off");
        names.push_back("D:/Codes/Boolean2016/exps/data/cmp_meshworks/lion.off");
        corkeval(names, expr, "D:/result.off");
        system("pause");
        return 0;
    }

    std::ifstream configFile(config_filename);
    if (!configFile.is_open())
    {
        std::cerr << "IO error: " << config_filename << std::endl;
        return -1;
    }
    int nMesh;
    std::string name, output, line;
    std::stringstream buffer;

    try
    {
        while (1)
        {
            nMesh = 0;
            names.clear();
            buffer = getValidStream(configFile);
            buffer >> nMesh;
            if (nMesh == 0 || configFile.eof()) break;

            for (int j = 0; j < nMesh; j++)
            {
                buffer = getValidStream(configFile);
                buffer >> name;
                if (name[0] == '@')
                {
                    --j;
                    name = name.substr(1, name.length());
                    buffer = getValidStream(configFile);
                    int s, e;
                    char nameBuffer[512];
                    buffer >> s >> e;
                    for (int k = s; k <= e; k++, j++)
                    {
                        sprintf(nameBuffer, name.c_str(), k);
                        names.push_back(nameBuffer);
                    }
                }
                else names.push_back(name);
            }
            // expr
            buffer = getValidStream(configFile);
            char first = buffer.peek();
            if (first == '@')
            {
                buffer.get();
                expr = "0";
                int num, count = 1;
                std::string opstr;
                char opch;
                while (1)
                {
                    buffer >> num >> opstr;
                    assert(opstr.size() == 1);
                    opch = opstr[0];
                    for (int j = 0; j < num; j++, count++)
                    {
                        expr += opch;
                        expr += std::to_string(count);
                    }
                    if (buffer.eof()) break;
                }
            }
            else buffer >> expr;

            // output
            buffer = getValidStream(configFile);
            buffer >> output;

            corkeval(names, expr, output);
        }
    }
    catch (...)
    {
        std::cerr << "non-standard ends or error occurred.\n";
    }

    configFile.close();

    if (argc == 1)
        system("pause");
    return 0;
}





