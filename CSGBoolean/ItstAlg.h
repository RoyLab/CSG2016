#pragma once

#include <unordered_map>
#include <unordered_set>
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

        typedef MyMesh::Face_handle FH;
        typedef std::unordered_set<IndexPair> TriIdSet;
        typedef std::unordered_map<IndexPair, TriIdSet*> MeshIdTriIdMap;
    public:
        ItstAlg(std::vector<MyMesh*>* meshes);
        ~ItstAlg();

        void doIntersection(std::vector<Octree::Node*>& intersectLeaves);
        bool IntersectionTest(FH fh0, FH fh1, TriIdSet* overlaps);
        bool checkManifoldEdge(FH fh0, FH fh1, TriIdSet* overlaps, TriTriIsectResult<K> &result, bool res[]);

    private:

        /*  合并优先级：
        已经登记为共享点的，按顺序排大小，先到大
        */
        void mergeProxy(VProxies::iterator a, VProxies::iterator b);

        VProxies                vProxy;
        VEntities               vEnt;

    };

}

