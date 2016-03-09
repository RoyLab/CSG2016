#pragma once

#include <vector>

#include "macroutil.h"
#include "AlgUserData.h"
#include "COctree.h"
#include "UndirectedGraph.h"

namespace CSG
{
    class AdjacentGraph :
        public myext::UndirectedGraph<bool>
    {
    public:
        void getIntersectPrimitives(int meshId, std::vector<int>& prims);
    };


    class ItstAlg
    {
        COMMON_PROPERTY_POINTER(AdjacentGraph, adjGraph);
        COMMON_PROPERTY_POINTER(std::vector<MyMesh*>, pMeshList);
    public:
        ItstAlg(std::vector<MyMesh*>* meshes);
        ~ItstAlg();

        void doIntersection(std::vector<Octree::Node*>& intersectLeaves);

    private:

        /*  合并优先级：
        已经登记为共享点的，按顺序排大小，先到大
        */
        void mergeProxy(VProxies::iterator a, VProxies::iterator b);

        VProxies                vProxy;
        VEntities               vEnt;

    };

}

