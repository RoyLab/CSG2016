#pragma once
#include <string>
#include <vector>
#include <queue>
#include <boost\shared_ptr.hpp>

#include "Octree.h"
#include "CGALext.h"
#include "MyMesh.h"

namespace CSG
{
    class ItstAlg;

    typedef K::Vector_3 Vec3d;
    typedef K::Point_3  Point3d;
    typedef K::Ray_3    Ray;

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
            boost::shared_ptr<IIndicatorVector> indicators;
            virtual ~SeedInfo() {}
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
            //boost::scoped_ptr<TrimCSGTree<MyMesh>> ttree1, ttree2;
        };

    public:
        MyAlgorithm(){}
        ~MyAlgorithm(){}

        void solve(const std::string& expr, std::vector<MyMesh*>& meshList);
        MyMesh* getResultMesh();
        MyMesh* popResultMesh();

    private:
        void floodColoring(CSGTree<MyMesh>* pCsg, ItstAlg* itstAlg);
        void createFirstSeed(SeedInfoWithId& info);
        IIndicatorVector* computeFullIndicator(VH fh, size_t meshId);
        void floodComplexGroup(GroupParseInfo& infos, SeedInfoWithHint& s);
        void floodSimpleGroup(GroupParseInfo& infos, SeedInfoWithHint& s);

    private:

        MyMesh*                     csgResult = nullptr;
        std::vector<MyMesh*>*       pMeshList = nullptr;

        Octree *                    pOctree = nullptr;
    };

    Relation PolyhedralInclusionTest(Point3d& point, Octree* pOctree, std::vector<MyMesh*>& pMesh, unsigned meshId, bool IsInverse);

}

