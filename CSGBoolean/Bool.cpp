#include <string>

#define CSG_EXPORTS

#include "bool.h"
#include "MyAlgorithm.h"

namespace CSG
{
    extern "C" CSG_API int loadMesh(std::vector<MyMesh*>& meshes, const std::vector<std::string>& names)
    {
        return 0;
    }


    void test()
    {
        std::vector<std::string> names;
        std::string expr;

        std::vector<MyMesh*> meshList;
        loadMesh(meshList, names);

        MyAlgorithm *alg = new MyAlgorithm;
        alg->solve(expr, meshList);

        MyMesh* res = alg->popResultMesh();
        SAFE_DELETE(alg);
    }
}