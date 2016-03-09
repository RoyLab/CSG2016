#pragma once
#include <string>
#include <vector>
#include <boost\smart_ptr.hpp>

#include "COctree.h"
#include "CGALext.h"
#include "MyMesh.h"
#include "TrimCSGTree.h"

namespace CSG
{
    class MyAlgorithm
    {
        typedef MyMesh::Face_handle     FH;
        typedef MyMesh::Vertex_handle   VH;

        template <class T>
        class Queue : public std::queue < T, std::list<T> > {};

        struct SeedInfo
        {
            FH              seedFacet;
            VH              seedVertex;
            IndicatorVector *indicators;
            size_t          meshId = -1;

            ~SeedInfo() { SAFE_DELETE(indicators); }
        };

        struct GroupParseInfo
        {
            Queue<SeedInfo> otherMeshSeedQueue;
            Queue<SeedInfo> curMeshSeedQueue;

            int32_t             curMeshId = -1;
            boost::scoped_ptr<TrimCSGTree<MyMesh>> ttree1, ttree2;
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
        void createFirstSeed(SeedInfo& info);
        IndicatorVector* computeFullIndicator(VH fh, size_t meshId);
        void floodComplexGroup();
        void floodSimpleGroup();

    private:

        MyMesh*                     csgResult = nullptr;
        std::vector<MyMesh*>*       pMeshList = nullptr;
    };

}

