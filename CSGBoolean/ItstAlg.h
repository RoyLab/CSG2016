#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "macroutil.h"
#include "AlgUserData.h"
#include "Octree.h"
#include "plane_reps.h"
#include "UndirectedGraph.h"

namespace CSG
{
    enum Sign
    {
        UNKOWN = 0,
        INTERSECT_ON_LINE,
        INTERSECT_ON_POINT,
        COPLANAR,
        LESS_THAN_INTERSECT_ON_POINT,
        NOT_INTERSECT
    };

    template <class _R>
    struct TriTriIsectResult
    {
        PosTag tagA[2], tagB[2];
        PBPoint<_R> A, B;
    };

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
        bool IntersectionTest(FH fh0, FH fh1, TriIdSet* overlaps, uint32_t meshId[2]);
        bool checkManifoldEdge(FH fh0, FH fh1, TriIdSet* overlaps, TriTriIsectResult<K> &result, int res[], uint32_t meshId[2]);

        /*  合并优先级：
        已经登记为共享点的，按顺序排大小，先到大
        */
        void mergeProxy(VProxyItr a, VProxyItr b);
        int checkDuplicatedPoints(PBPoint<K>& point, FH fhs, PosTag tags, VProxyItr& outcome);
        int checkDuplicatedPoints(std::vector<VProxyItr>& plist, PBPoint<K>& point, VProxyItr& proxy);
        void getVProxy(PBPoint<K>& point, int addwhat[2], FH fhs[2], PosTag tags[2], int oId[2], VProxyItr outproxy[2], uint32_t meshId[2]);
        int addPoint(FH fh, PosTag tags, VProxyItr proxy);
        VProxyItr addVEntity(PBPoint<K>& point);
        void resolveIntersection(FH tri, size_t meshId);

        VProxies                vProxy;
        VEntities               vEnt;
    };

}

