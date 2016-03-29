
#include "BSPTree.h"
#include "plane_reps.h"

namespace CSG
{
    BSPTree::BSPTree(std::vector<FH>& facets)
    {
        std::vector<Polygon> polygons;
        for (auto &face : facets)
        {
            if (!face->data->planeRep)
            {
                face->data->planeRep = new PBTriangle<K>(face->data->triangle);
                polygons.emplace_back(*face->data->planeRep);
            }
        }

        buildBSP(polygons);
    }

    BSPTree::Node* BSPTree::buildBSP(std::vector<Polygon>& facets)
    {
        Node* pNode = new Node;
        pNode->sp = facets[0].get_sp();

        std::vector<Polygon> fronts;
        std::vector<Polygon> backs;

        splitPloygonsWithPlane(facets, pNode->sp, fronts, backs, pNode->coins);

        auto dir1 = pNode->sp.orthogonal_direction(),
            dir2 = pNode->coins[0].get_sp().orthogonal_direction();
        bool bSameOrient = dir1.to_vector() * dir2.to_vector() > 0.0;
        if (fronts.empty())
        {
            pNode->left = new Node;
            if (bSameOrient)
                pNode->left->relation = In;
            else
                pNode->left->relation = Out;
        }
        else
            pNode->left = buildBSP(fronts);
        if (backs.empty())
        {
            pNode->right = new Node;
            if (bSameOrient)
                pNode->right->relation = Out;
            else
                pNode->right->relation = In;
        }
        else
            pNode->right = buildBSP(backs);
        return pNode;
    }

    void BSPTree::splitPloygonsWithPlane(std::vector<Polygon>& facets, Plane& sp,
        std::vector<Polygon>& fronts, std::vector<Polygon>& backs, std::vector<Polygon>& coins)
    {
        Polygon front, back;
        fronts.clear();
        backs.clear();
        coins.clear();

        for (int i = 0; i < facets.size(); i++)
        {
            front.Clear();
            back.Clear();
            RelationToPlane rp = facets[i].ClipByPlane(sp, front, back);
            if (rp == On)
                coins.push_back(facets[i]);
            else if (rp == Straddling)
            {
                fronts.push_back(front);
                backs.push_back(back);
            }
            else if (rp == Front)
                fronts.push_back(front);
            else
                backs.push_back(back);
        }
    }

}