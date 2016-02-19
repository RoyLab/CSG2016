#pragma once
#include <string>
#include <vector>
#include "COctree.h"
#include "CGALext.h"

namespace CSG
{
    class MyMesh;

    typedef MyMesh::Vertex_iterator     PointListItr;
    typedef std::list<PointListItr>     PointListItrList;
    typedef PointListItrList::iterator  PointListItrListItr;

    struct SharedEdge
    {
        std::vector<PointListItrListItr>  innerPoints;
    };

    struct IsectTriangleInfo
    {
        std::vector<PointListItrListItr>  innerPoints;
        std::vector<MyMesh::Face_handle> coplanars;
    };


    class MyAlgorithm
    {
    public:
        MyAlgorithm(){}
        ~MyAlgorithm(){}

        void solve(const std::string& expr, std::vector<MyMesh*>& meshList);
        MyMesh* getResultMesh();
        MyMesh* popResultMesh();

    private:
        void doIntersection(std::vector<MyMesh*>& meshList, std::vector<Octree::Node*>& intersectLeaves);
        void floodColoring(CSGTree<MyMesh>* pCsg, Octree* pOctree);
        void checkNonmanifoldEdge(MyMesh::Face_handle, MyMesh::Face_handle, myext::TriTriIsectResult<K>*, void*);

        void setupIsectFacet(MyMesh::Face_handle fh);

        /*  合并优先级：
        已经登记为共享点的，按顺序排大小，先到大
        */
        void setupPonits(MyMesh::Face_handle fh0, MyMesh::Face_handle fh1, const myext::TriTriIsectResult<K>& result);

    private:
        MyMesh*                     csgResult = nullptr;
        PointListItrList            pointAgency;
        //SharedZone*                 szone = nullptr;
    };

}

