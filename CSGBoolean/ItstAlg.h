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

        /*  �ϲ����ȼ���
        �Ѿ��Ǽ�Ϊ�����ģ���˳���Ŵ�С���ȵ���
        */
        void mergeProxy(VProxies::iterator a, VProxies::iterator b);

        VProxies                vProxy;
        VEntities               vEnt;

    };

}

