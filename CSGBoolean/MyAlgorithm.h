#pragma once
#include <string>
#include <vector>
#include "COctree.h"

namespace CSG
{
    class MyMesh;

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

    private:
        MyMesh* result = nullptr;
    };

}

