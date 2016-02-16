#include "MyAlgorithm.h"
#include "csg.h"
#include "COctree.h"
#include "MPMesh.h"


namespace CSG
{
    void MyAlgorithm::solve(const std::string& expr, std::vector<MPMesh*>& meshes)
    {
        std::vector<MPMesh*> meshList;
        for (MPMesh* pMesh : meshes)
            meshList.push_back(new MPMesh(*pMesh));

        CSGTree<MPMesh>* pCsg = new CSGTree<MPMesh>;
        pCsg->createCSGTreeFromExpr(expr, meshList.data(), meshList.size());
        pCsg->makePositiveAndLeftHeavy();

        Octree *pOctree = new Octree;
        pOctree->build(meshList, &octreeLeaves);

        doIntersection();
        floodColoring(pCsg, pOctree);

        SAFE_DELETE(pCsg);
        SAFE_DELETE(pOctree);
        for (MPMesh* pMesh : meshList)
            delete pMesh;
    }
}