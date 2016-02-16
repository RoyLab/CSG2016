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
        void doIntersection();
        void floodColoring(CSGTree<MyMesh>* pCsg, Octree* pOctree);

    private:
        std::vector<Octree::Node*> octreeLeaves;
        MyMesh* result = nullptr;
    };

}

