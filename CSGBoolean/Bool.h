#pragma once
#include <vector>


#ifdef CSG_EXPORTS
#define CSG_API __declspec(dllexport)
#else
#define CSG_API __declspec(dllimport)
#endif

namespace CSG
{
    class MyMesh;

    extern "C" CSG_API int loadMesh(std::vector<MyMesh*>& meshes, const std::vector<std::string>& names);
}

