#pragma once

#include <vector>

#include "AlgUserData.h"
#include "COctree.h"
#include "UndirectedGraph.h"

namespace CSG
{

    class ItstAlg
    {
        typedef myext::UndirectedGraph<bool> AdjacentGraph;

    public:
        ItstAlg(std::vector<MyMesh*>* meshes);
        ~ItstAlg();

        void doIntersection(std::vector<Octree::Node*>& intersectLeaves);

    private:

        /*  合并优先级：
        已经登记为共享点的，按顺序排大小，先到大
        */
        void mergeProxy(VProxies::iterator a, VProxies::iterator b);

        // ref
        std::vector<MyMesh*>    *pMeshList = nullptr;

        //entity
        AdjacentGraph           *adjGraph = nullptr;

        VProxies                vProxy;
        VEntities               vEnt;

    };

}

