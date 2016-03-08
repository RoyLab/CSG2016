#pragma once
#include <string>
#include <vector>
#include <boost\smart_ptr.hpp>

#include "COctree.h"
#include "CGALext.h"
#include "MyMesh.h"

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
            Queue<SeedInfoWithMeshId> otherMeshSeedQueue;
            Queue<SeedInfo> curMeshSeedQueue;
            int32_t curMeshId = -1;
        };

    public:
        MyAlgorithm(){}
        ~MyAlgorithm(){}

        void solve(const std::string& expr, std::vector<MyMesh*>& meshList);
        MyMesh* getResultMesh();
        MyMesh* popResultMesh();

    private:
        void floodColoring(CSGTree<MyMesh>* pCsg, ItstAlg* itstAlg);
        void copyAutoIndicator(AutoIndicator& target, AutoIndicator& source);
        void createFirstSeed(SeedInfoWithMeshId& info);
        AutoIndicator computeFullIndicator(VH fh, size_t meshId);
        void floodComplexGroup();
        void floodSimpleGroup();

    private:

        MyMesh*                     csgResult = nullptr;
        std::vector<MyMesh*>*       pMeshList = nullptr;
    };

}

