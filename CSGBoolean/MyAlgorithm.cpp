#include "MyAlgorithm.h"
#include "csg.h"
#include "COctree.h"
#include "MyMesh.h"


namespace CSG
{
    void MyAlgorithm::solve(const std::string& expr, std::vector<MyMesh*>& meshes)
    {
        std::vector<MyMesh*> meshList;
        for (MyMesh* pMesh : meshes)
            meshList.push_back(new MyMesh(*pMesh));

        CSGTree<MyMesh>* pCsg = new CSGTree<MyMesh>;
        pCsg->createCSGTreeFromExpr(expr, meshList.data(), meshList.size());
        pCsg->makePositiveAndLeftHeavy();

        Octree *pOctree = new Octree;
        pOctree->build(meshList, &octreeLeaves);

        doIntersection();
        floodColoring(pCsg, pOctree);

        SAFE_DELETE(pCsg);
        SAFE_DELETE(pOctree);
        for (MyMesh* pMesh : meshList)
            delete pMesh;
    }
}