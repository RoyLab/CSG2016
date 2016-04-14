
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
            }
            polygons.emplace_back(*face->data->planeRep, face);
        }

        mp_root = buildBSP(polygons);
    }
    
    Relation BSPTree::determine(const PBPoint<K>& point)
    {
        return determine(point, mp_root);
    }

    Relation BSPTree::determine(const PBPoint<K>& point, Node* node)
    {
        if (node->relation == PR_In)
            return REL_INSIDE;

        if (node->relation == PR_Out)
            return REL_OUTSIDE;

        assert(node->relation == PR_None);
        assert(node->left && node->right);

        RelationToPlane rel = point.classifyByPlane(node->sp); 
        switch (rel)
        {
        case CSG::Front:
            return determine(point, node->left);
        case CSG::On:
        {
            if (node->coins.size() > 1) ReportError("more than 1 coin, not strange.");
            Relation rell = determine(point, node->left);
            Relation relr = determine(point, node->right);
            if (rell == relr) return rell;
            else
            {
                Plane_ext<K> k;
                for (auto& polys : node->coins)
                {
                    assert(polys.fh->data->planeRep);
                    if (polys.fh->data->planeRep->insideBPs(point))
                        m_lastCoins.insert(node->coins[0].fh);
                }
                return REL_ON_BOUNDARY;
            }
        }
        case CSG::Behind:
            return determine(point, node->left);
        default:
            assert(0);
            return REL_NOT_AVAILABLE;
        }
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
                pNode->left->relation = PR_Out;
            else
                pNode->left->relation = PR_In;
        }
        else
            pNode->left = buildBSP(fronts);
        if (backs.empty())
        {
            pNode->right = new Node;
            if (bSameOrient)
                pNode->right->relation = PR_In;
            else
                pNode->right->relation = PR_Out;
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