#include "precompile.h"
#include "hybrid_geometry.h"

#include "RegularMesh.h"
#include "BSP.h"

namespace Boolean
{
    void BSPTree::buildNoCross(std::vector<Triangle*>& faces)
    {
        mp_root = buildRecursion(faces);
    }

    void RecursiveRelease(BSPTree::Node* node)
    {
        if (!node->isLeaf())
        {
            RecursiveRelease(node->negChild);
            RecursiveRelease(node->posChind);
        }
        delete node;
    }

    void BSPTree::reset()
    {
        RecursiveRelease(mp_root);
        mp_root = nullptr;
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

    BSPTree::Node * BSPTree::buildRecursion(std::vector<Triangle*>& faces)
    {
        Node* pNode = new Node;
        pNode->sp = faces[0]->supportingPlane();

        std::vector<Triangle*> fronts;
        std::vector<Triangle*> backs;

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