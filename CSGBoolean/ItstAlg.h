#pragma once

#include <vector>

#include "macroutil.h"
#include "AlgUserData.h"
#include "Octree.h"
#include "UndirectedGraph.h"

namespace CSG
{
    class AdjacentGraph :
        public myext::UndirectedGraph<bool>
    {
    public:
        AdjacentGraph(size_t n) :myext::UndirectedGraph<bool>(n){}
        void getIntersectPrimitives(int meshId, std::vector<int>& prims);
    };


    class ItstAlg
    {
        COMMON_PROPERTY_POINTER(AdjacentGraph, adjGraph);
        COMMON_PROPERTY_POINTER(std::vector<MyMesh*>, meshList);
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

