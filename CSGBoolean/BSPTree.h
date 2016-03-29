#pragma once

#include <vector>

#include "MyMesh.h"
#include "plane_reps.h"

namespace CSG
{
    class BSPTree
    {
        typedef Plane_ext<K> Plane;
        typedef MyMesh::Face_handle FH;
        typedef PBPolygon<K> Polygon;

    private:
        struct Node
        {
            Plane sp;
            std::vector<Polygon> coins;
            std::vector<Plane> bounds;
            PosRelation   relation;
            Node *left = nullptr, *right = nullptr;
        };

    public:
        BSPTree(std::vector<FH>& facets);

    private:
        Node* buildBSP(std::vector<Polygon>& facets);
        void splitPloygonsWithPlane(std::vector<Polygon>& facets, Plane& sp,
            std::vector<Polygon>& fronts, std::vector<Polygon>& backs, std::vector<Polygon>& coins);

        Node* mp_root = nullptr;

    };
}