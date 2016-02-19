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

        /*  �ϲ����ȼ���
        �Ѿ��Ǽ�Ϊ�����ģ���˳���Ŵ�С���ȵ���
        */
        void setupPonits(MyMesh::Face_handle fh0, MyMesh::Face_handle fh1, const myext::TriTriIsectResult<K>& result);
        void copyAutoIndicator(AutoIndicator& target, AutoIndicator& source);
        void createFirstSeed(SeedInfoWithMeshId& info);
        AutoIndicator computeFullIndicator(VH fh);


    private:
        MyMesh*                     csgResult = nullptr;
        std::vector<MyMesh*>*       pMeshList = nullptr;

        myext::TriangleTable<bool>* meshRelTable = nullptr;
        PointListItrList            pointAgency;
    };

}

