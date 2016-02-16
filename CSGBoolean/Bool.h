#pragma once
#include <vector>


#ifdef CSG_EXPORTS
#define CSG_API __declspec(dllexport)
#else
#define CSG_API __declspec(dllimport)
#endif

namespace CSG
{
 //   struct MPMesh2;
 //   extern "C" CSG_API GS::BaseMesh* BooleanOperation(GS::CSGExprNode* input, void* stdoutput, bool db);
    //extern "C" CSG_API GS::BaseMesh* BooleanOperation_MultiThread(GS::CSGExprNode* input);
 //   extern "C" CSG_API MPMesh2* BooleanOperationFeito(MPMesh2* mesh1,  MPMesh2* mesh2, int op, void* hd);
 //   extern "C" CSG_API MPMesh2* Convert2MPMesh2(const GS::BaseMesh* mesh);
 //   extern "C" CSG_API GS::BaseMesh* Convert2BaseMesh2(MPMesh2* mesh);
    //extern "C" CSG_API void SnapModel(const GS::BaseMesh* mesh);
 //   extern "C" CSG_API void ReleaseMPMesh2(MPMesh2*& ptr);

    class MyMesh;

    extern "C" CSG_API int loadMesh(std::vector<MyMesh*>& meshes, const std::vector<std::string>& names);
}

