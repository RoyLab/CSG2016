#pragma once
#include <string>
#include <vector>
#include "COctree.h"

namespace CSG
{
    class MPMesh;

    class MyAlgorithm
    {
    public:
        MyAlgorithm(){}
        ~MyAlgorithm(){}

        void solve(const std::string& expr, std::vector<MPMesh*>& meshList);
        MPMesh* getResultMesh();
        MPMesh* popResultMesh();

    private:
        void doIntersection();
        void floodColoring(CSGTree<MPMesh>* pCsg, Octree* pOctree);

    private:
        std::vector<Octree::Node*> octreeLeaves;
        MPMesh* result = nullptr;
    };

}

