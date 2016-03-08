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

        /*  �ϲ����ȼ���
        �Ѿ��Ǽ�Ϊ�����ģ���˳���Ŵ�С���ȵ���
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

