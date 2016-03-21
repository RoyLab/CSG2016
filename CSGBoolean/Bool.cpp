#include <string>
#include <fstream>
#include <CGAL/IO/Polyhedron_iostream.h>

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
            meshes.push_back(newPoly);
        }

        return 0;
    }


    extern "C" CSG_API void test()
    {
        GS::exactinit();
        std::vector<std::string> names;
        names.push_back("../../models/box1.off");
        names.push_back("../../models/box2.off");

        std::string expr("0+1");

        std::vector<MyMesh*> meshList;
        loadMesh(meshList, names);

        MyAlgorithm *alg = new MyAlgorithm;
        alg->solve(expr, meshList);

        for (auto mesh : meshList)
        {
            //std::cout << mesh->size_of_facets() << std::endl;
            SAFE_DELETE(mesh);
        }

        SAFE_DELETE(alg);
    }
}