#pragma once
#include <string>
#include <vector>
#include <boost\shared_ptr.hpp>

#include "COctree.h"
#include "CGALext.h"
#include "MyMesh.h"
#include "TriangleTable.h"

namespace CSG
{
    class MyAlgorithm
    {
        typedef MyMesh::Face_handle     FH;
        typedef MyMesh::Vertex_handle   VH;
        typedef int8_t                  Indicator;
        typedef boost::shared_ptr<Indicator[]> AutoIndicator;

        template <class T>
        class Queue : public std::queue < T, std::list<T> > {};

        struct SeedInfo
        {
            FH              seedFacet;
            VH              seedVertex;
            AutoIndicator   indicator;
        };

        struct SeedInfoWithMeshId :
            public SeedInfo
        {
            size_t          meshId;
        };

        struct GroupParseInfo
        {
            Queue<SeedInfoWithMeshId > otherMeshSeedQueue;
            Queue<SeedInfo > curMeshSeedQueue;
            size_t curMeshId = -1;
            bool *meshSeedFlag = nullptr;
        };

    public:
        MyAlgorithm(){}
        ~MyAlgorithm(){}

        void solve(const std::string& expr, std::vector<MyMesh*>& meshList);
        MyMesh* getResultMesh();
        MyMesh* popResultMesh();

    private:
        void doIntersection(std::vector<Octree::Node*>& intersectLeaves);
        void floodColoring(CSGTree<MyMesh>* pCsg);
        void checkNonmanifoldEdge(MyMesh::Face_handle, MyMesh::Face_handle, myext::TriTriIsectResult<K>*, void*);

        void setupIsectFacet(MyMesh::Face_handle fh);

        /*  合并优先级：
              已经登记为共享点的，按顺序排大小，先到大
        */
        void setupPonits(MyMesh::Face_handle fh0, MyMesh::Face_handle fh1, const myext::TriTriIsectResult<K>& result);
        void copyAutoIndicator(AutoIndicator& target, AutoIndicator& source);
        void createFirstSeed(SeedInfoWithMeshId& info);
        AutoIndicator computeFullIndicator(VH fh, size_t meshId);

        void floodComplexGroup();
        void floodSimpleGroup();

    private:

        MyMesh*                     csgResult = nullptr;
        std::vector<MyMesh*>*       pMeshList = nullptr;

        myext::TriangleTable<bool>* meshRelTable = nullptr;
        PointListItrList            pointAgency;
    };

}

