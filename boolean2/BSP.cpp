#include "precompile.h"
#include "RegularMesh.h"
#include "BSP.h"

namespace Boolean
{
    Oriented_side orientation(const XPlane& p, const MyVertex& v);

    void BSPTree::buildNoCross(std::vector<IPolygon*>& faces)
    {
        mp_root = buildRecursion(faces);
    }

    void BSPTree::build(std::vector<IPolygon*>& faces)
    {
        assert(0);
    }

    void BSPTree::reset()
    {
        SAFE_DELETE(mp_root);
    }

    Oriented_side BSPTree::classify(const MyVertex & v, XPlane * bspPlane) const
    {
        assert(mp_root);

        Oriented_side side;
        Node* curNode = mp_root;
        while (!curNode->isLeaf())
        {
            side = orientation(curNode->sp, v);

            if (side == ON_POSITIVE_SIDE)
                curNode = curNode->posChind;
            else if (side == ON_NEGATIVE_SIDE)
                curNode = curNode->negChild;
            else
            {
                if (bspPlane)
                    *bspPlane = curNode->sp;
                return ON_ORIENTED_BOUNDARY;
            }
        }
        return curNode->relation;
    }

    BSPTree::Node * BSPTree::buildRecursion(std::vector<IPolygon*>& faces)
    {
        Node* pNode = new Node;
        pNode->sp = faces[0]->supportingPlane();

        std::vector<IPolygon*> fronts;
        std::vector<IPolygon*> backs;

        size_t i, j; Oriented_side side;
        for (i = 0; i < faces.size(); i++)
        {
            for (j = 0; j < faces[i]->degree(); j++)
            {
                side = orientation(pNode->sp, faces[i]->vertex(j));
                if (side == ON_POSITIVE_SIDE)
                {
                    fronts.push_back(faces[i]);
                    break;
                }
                
                if (side == ON_NEGATIVE_SIDE)
                {
                    backs.push_back(faces[i]);
                    break;
                }
                //else pNode->coins.push_back(faces[i]);
            }
        }

        if (fronts.empty())
        {
            pNode->posChind = new Node;
            pNode->posChind->relation = ON_POSITIVE_SIDE;
        }
        else 
            pNode->posChind = buildRecursion(fronts);

        if (backs.empty())
        {
            pNode->negChild = new Node;
            pNode->negChild->relation = ON_NEGATIVE_SIDE;
        }
        else
            pNode->negChild = buildRecursion(backs);

        return pNode;
    }


}