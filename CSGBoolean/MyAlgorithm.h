#pragma once
#include <string>
#include <vector>
#include <queue>
#include <boost\smart_ptr.hpp>

#include "Octree.h"
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
        class Queue: public std::queue<T,std::list<T>> {};

        struct SimpleSeedInfo
        {
            FH  seedFacet;
            VH  seedVertex;
        };

        struct SeedInfo:
            public SimpleSeedInfo
        {
            IndicatorVector *indicators;
            virtual ~SeedInfo() { SAFE_DELETE(indicators); }
        };

        struct SeedInfoWithId:
            public SeedInfo
        {
            size_t          meshId = -1;
        };

        class Hint
        {
            COMMON_PROPERTY(uint8_t, source);
            COMMON_PROPERTY(uint8_t, target);
        public:
            Hint() { makeInvalid(); }
            Hint(uint8_t s, uint8_t t): m_source(s), m_target(t) {}
            bool isValid() const { return m_source == m_target; }
            void makeInvalid() { m_source = 0; m_target = 0; }
        };

        struct SeedInfoWithHint :
            public SeedInfo
        {
            Hint hint;
        };

        struct SimpleSeedInfoWithHint :
            public SimpleSeedInfo
        {
            Hint hint;
        };

        struct GroupParseInfo
        {
            Queue<SeedInfoWithId> otherMeshSeedQueue;
            Queue<SeedInfoWithHint> curMeshSeedQueue;

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
        void createFirstSeed(SeedInfoWithId& info);
        IndicatorVector* computeFullIndicator(VH fh, size_t meshId);
        void floodComplexGroup(GroupParseInfo& infos, SeedInfoWithHint& s);
        void floodSimpleGroup(GroupParseInfo& infos, SeedInfoWithHint& s);

    private:

        MyMesh*                     csgResult = nullptr;
        std::vector<MyMesh*>*       pMeshList = nullptr;

        Octree *                    pOctree = nullptr;
    };

    Relation pointInPolyhedron(CGAL::Point_3<K>& p, MyMesh* mesh, Octree* pOctree);

}

