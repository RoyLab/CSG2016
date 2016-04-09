#pragma once
#include <string>
#include <vector>
#include <queue>
#include <boost\shared_ptr.hpp>

#include "Octree.h"
#include "CGALext.h"
#include "MyMesh.h"
#include "AlgUserData.h"
#include "csg.h"
#include "ItstGraph.h"

namespace CSG
{
    class ItstAlg;

    typedef K::Vector_3 Vec3d;
    typedef K::Point_3  Point3d;
    typedef K::Ray_3    Ray;

    template <class HDS> class Delegate;

    class MyAlgorithm
    {
        typedef MyMesh::Face_handle     FH;
        typedef MyMesh::Vertex_handle   VH;
        typedef MyMesh::Halfedge_handle   EH;
        typedef MyMesh::HalfedgeDS  HalfedgeDS;

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
            Queue<SeedInfo> curMeshSeedQueue;

            int32_t             curMeshId = -1;
            boost::shared_ptr<CSGTreeNode> ttree1, ttree2;
            CSGTreeNode** curTreeLeaves = nullptr;
        };

    public:
        MyAlgorithm(){}
        ~MyAlgorithm(){}

        void solve(const std::string& expr, std::vector<MyMesh*>& meshList);
        boost::shared_ptr<MyMesh> getResultMesh() { return csgResult; }

    private:
        void floodColoring(CSGTree<MyMesh>* pCsg, ItstAlg* itstAlg);
        void createFirstSeed(SeedInfoWithId& info);
        IIndicatorVector* computeFullIndicator(VH fh, size_t meshId);
        void floodComplexGroup(GroupParseInfo& infos, SeedInfo& s);
        void floodSimpleGroup(GroupParseInfo& infos, SeedInfo& s);
        void figureOutFaceInds(SeedInfo& s, int meshId);
        void figureOutFaceInds(VH p, VH q, VH r, int meshId, IIndicatorVector* inds);

        void AddFacet(FH fh);
        bool isSameGroup(FH fh0, FH fh1) const;
        void genVertexInds(IIndicatorVector* target, VH vh) const;
        bool needAdd(FH fh, ItstGraph::Loop& loop, TestTree& testList);
        void addLoop(ItstGraph::Loop&);

    private:
        boost::shared_ptr<MyMesh>   csgResult;
        std::vector<MyMesh*>        *pMeshList = nullptr;
        Delegate<HalfedgeDS>        *pConstruct = nullptr;
        Octree                      *pOctree = nullptr;
        ItstAlg                     *itst = nullptr;
    };

    Relation PolyhedralInclusionTest(Point3d& point, Octree* pOctree, std::vector<MyMesh*>& pMesh, unsigned meshId, bool IsInverse);

    template <class HDS>
    class Delegate : public CGAL::Modifier_base<HDS> {
        
        struct Tripple
        { 
            Tripple(int i, int j, int k){ idx[0] = i; idx[1] = j; idx[2] = k; }
            int idx[3]; 
        };

    public:
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;

        Delegate() {}

        void operator()(HDS& hds) {

            CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);
            B.begin_surface(points.size(), indices.size() / 3);

            for (size_t i = 0; i < points.size(); i++)
                B.add_vertex(points[i]);

            for (size_t i = 0; i < indices.size(); i++)
            {
                B.begin_facet();
                for (int j = 0; j < 3; j++)
                    B.add_vertex_to_facet(indices[i].idx[j]);
                B.end_facet();
            }

            B.end_surface();
        }

        size_t add_vertex(const Point& pts)
        {
            points.push_back(pts);
            return points.size() - 1;
        }

        void addFacets(int i, int j, int k)
        {
            indices.emplace_back(i, j, k);
        }

    private:
        std::deque<Point> points;
        std::deque<Tripple> indices;
    };
}

