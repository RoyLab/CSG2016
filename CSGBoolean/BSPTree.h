#pragma once
#include <set>
#include <vector>

#include "csgdefs.h"
#include "MyMesh.h"
#include "plane_reps.h"

namespace CSG
{
    class BSPTree
    {
        typedef Plane_ext<K> Plane;
        typedef MyMesh::Face_handle FH;
        struct Polygon :
            public PBPolygon<K>
        {
            Polygon(){}
            Polygon(const PBTriangle<K>& triangle, FH f) :
                PBPolygon<K>(triangle), fh(f){}
            FH fh;

            RelationToPlane ClipByPlane(const PlaneExt& bp, Polygon& front, Polygon& back) const
            {
                front.fh = fh;
                back.fh = fh;
                return PBPolygon<K>::ClipByPlane(bp, front, back);
            }
        };

    private:
        struct Node
        {
            Plane sp;
            std::vector<Polygon> coins;
            std::vector<Plane> bounds;
            PosRelation relation = PR_None;

            /* left·Åfronts */
            Node *left = nullptr, *right = nullptr;
        };

    public:
        BSPTree(std::vector<FH>& facets);
        Relation determine(const PBPoint<K>& point);
        const std::set<FH>& getCoins() const { return m_lastCoins; }

    private:
        Node* buildBSP(std::vector<Polygon>& facets);
        Relation determine(const PBPoint<K>& point, Node* node);
        void splitPloygonsWithPlane(std::vector<Polygon>& facets, Plane& sp,
            std::vector<Polygon>& fronts, std::vector<Polygon>& backs, std::vector<Polygon>& coins);

        Node* mp_root = nullptr;
        std::set<FH> m_lastCoins;
    };
}