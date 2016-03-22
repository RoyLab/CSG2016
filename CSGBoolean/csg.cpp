﻿#include "csg.h"
#include <stack>
//#include <cassert>
//#include <climits>
//#include <algorithm>

namespace CSG
{
    namespace
    {
        int OperatorPriority(int c)
        {
            if (c == ')')
                return 0;
            if ((c == '+') || (c == '-'))
                return 1;
            if (c == '*')
                return 2;
            return 3;// '('
        }
    }


    std::string InfixToPostfix(const std::string& infix)
    {
        std::string postfix;
        std::stack<char> opstack;
        for (int i = 0; i< infix.size(); i++)
        {
            if (!IsOperator(infix[i]))
            {
                postfix += infix[i];
                continue;
            }
            postfix += ' ';
            if (opstack.empty() || (opstack.top() == '(' && infix[i] != ')'))
                opstack.push(infix[i]);
            else {
                while (OperatorPriority(infix[i]) <= OperatorPriority(opstack.top()))
                {
                    if (opstack.top() != '(')
                    {
                        postfix += opstack.top();
                        opstack.pop();
                        if (opstack.empty())
                            break;
                    }
                    else
                        break;
                }
                if (infix[i] != ')')
                    opstack.push(infix[i]);
                else if (opstack.top() == '(')
                    opstack.pop();
            }
        }
        while (!opstack.empty())
        {
            postfix += opstack.top();
            opstack.pop();
        }
        return postfix;
    }



    int test()
    {
        //std::string str;
        //int **a = nullptr;

        //CSGTree<int> tree;
        //tree.createCSGTreeFromExpr(str, a, 3);
        return 0;
    }
    //unsigned long long mark;

//    static void ConvertCSGTreeNode(GS::CSGExprNode* input, CSGTree::Node** pRoot, std::vector<MPMesh*>& meshList)
//    {
//        if (!input) return;
//
//        CSGTree::Node*& root = *pRoot;
//        root = new CSGTree::Node;
//
//        GS::BaseMesh* pBaseMesh = input->GetMesh();
//        if (pBaseMesh)
//        {
//            root->pMesh = ConvertToMPMesh(pBaseMesh);
//            root->pMesh->ID = meshList.size();
//            root->Type = TYPE_LEAF;
//            meshList.emplace_back(root->pMesh);
//            return;
//        }
//
//        switch (input->BoolOperation())
//        {
//        case GS::CSGUnion:
//            root->Type = TYPE_UNION;
//            break;
//        case GS::CSGIntersect:
//            root->Type = TYPE_INTERSECT;
//            break;
//        case GS::CSGDiff:
//            root->Type = TYPE_DIFF;
//            break;
//        default:
//            root->Type = TYPE_UNKNOWN;
//            break;
//        }
//
//        ConvertCSGTreeNode(input->LeftNode(), &root->pLeft, meshList);
//        ConvertCSGTreeNode(input->RightNode(), &root->pRight, meshList);
//
//        if (root->pLeft) root->pLeft->Parent = root;
//        if (root->pRight) root->pRight->Parent = root;
//    }
//
//    CSGTree* ConvertCSGTree(GS::CSGExprNode* root, MPMesh*** arrMesh, int *nMesh) // convert nodes.
//    {
//        if (!root) return NULL;
//
//        CSGTree* pRes = new CSGTree;
//        pRes->pRoot = NULL;
//        std::vector<MPMesh*> meshList;
//
//        ConvertCSGTreeNode(root, &pRes->pRoot, meshList);
//
//        MPMesh**& pMesh = *arrMesh;
//        *nMesh = meshList.size();
//        pMesh = new MPMesh*[*nMesh];
//        memcpy(pMesh, meshList.data(), sizeof(MPMesh*)*(*nMesh));
//        return pRes;
//    }
//
//    static CSGTree::Node* ConvertToPositiveTree(const CSGTree::Node* root, bool inverse, unsigned level, unsigned& maxLvl)
//    {
//        CSGTree::Node* res = new CSGTree::Node;
//        if (root->Type == TYPE_LEAF) // it is a leaf node
//        {
//            res->bInverse = inverse;
//            res->Type = TYPE_LEAF;
//            res->pMesh = root->pMesh;
//            res->pMesh->bInverse = res->bInverse;
//            maxLvl = (maxLvl < level) ? level : maxLvl;
//        }
//        else
//        {
//            unsigned Ldepth(0), Rdepth(0);
//            if (root->Type == TYPE_DIFF)
//            {
//                res->Type = TYPE_INTERSECT;
//                res->pLeft = ConvertToPositiveTree(root->pLeft, inverse, level + 1, Ldepth);
//                res->pRight = ConvertToPositiveTree(root->pRight, !inverse, level + 1, Rdepth);
//            }
//            else
//            {
//                res->Type = root->Type;
//                res->pLeft = ConvertToPositiveTree(root->pLeft, inverse, level + 1, Ldepth);
//                res->pRight = ConvertToPositiveTree(root->pRight, inverse, level + 1, Rdepth);
//            }
//            maxLvl = (Ldepth > Rdepth) ? Ldepth : Rdepth;
//            if (Ldepth < Rdepth)
//            {
//                auto tmpNode = res->pLeft;
//                res->pLeft = res->pRight;
//                res->pRight = tmpNode;
//            }
//
//            if (res->pLeft) res->pLeft->Parent = res;
//            if (res->pRight) res->pRight->Parent = res;
//
//            if (inverse)
//                res->Type = (res->Type == TYPE_INTERSECT) ? TYPE_UNION : TYPE_INTERSECT;
//        }
//
//        return res;
//    }
//
//    void GetLeafList(CSGTree::Node* root, std::vector<int>& list)
//    {
//        if (root->Type == TYPE_LEAF)
//        {
//            list.push_back(root->pMesh->ID);
//        }
//        else
//        {
//            GetLeafList(root->pLeft, list);
//            GetLeafList(root->pRight, list);
//        }
//    }
//
//    static void GetLeafList(CSGTree::Node* root, decltype(CSGTree::Leaves)& leaves)
//    {
//        if (root->Type == TYPE_LEAF)
//        {
//            leaves[root->pMesh->ID] = root;
//        }
//        else
//        {
//            GetLeafList(root->pLeft, leaves);
//            GetLeafList(root->pRight, leaves);
//        }
//    }
//
//    static void GetLeafList(CSGTree* tree)
//    {
//        tree->Leaves.clear();
//        if (!tree->pRoot) return;
//        GetLeafList(tree->pRoot, tree->Leaves);
//    }
//
//    CSGTree* ConvertToPositiveTree(const CSGTree* myTree)
//    {
//        CSGTree* result = new CSGTree;
//        unsigned maxLvl = 0;
//        result->pRoot = ConvertToPositiveTree(myTree->pRoot, false, 0, maxLvl);
//        GetLeafList(result);
//        return result;
//    }
//
//    static CSGTree::Node* copy(const CSGTree::Node* thiz)
//    {
//        if (!thiz) return nullptr;
//
//        CSGTree::Node* pRes = new CSGTree::Node(*thiz);
//        pRes->pLeft = copy(thiz->pLeft);
//        pRes->pRight = copy(thiz->pRight);
//
//        if (pRes->pLeft) pRes->pLeft->Parent = pRes;
//        if (pRes->pRight) pRes->pRight->Parent = pRes;
//
//        return pRes;
//    }
//
//    CSGTree* copy(const CSGTree* thiz)
//    {
//        if (!thiz) return nullptr;
//
//        CSGTree* pCopy = new CSGTree;
//        pCopy->pRoot = copy(thiz->pRoot);
//        GetLeafList(pCopy);
//        return pCopy;
//    }
//
//    CSGTree::Node* copy2(const CSGTree::Node* thiz, CSGTree::Node** leafList)
//    {
//        if (!thiz) return nullptr;
//
//        CSGTree::Node* pRes = new CSGTree::Node(*thiz);
//
//        if (thiz->pMesh) leafList[thiz->pMesh->ID] = pRes;
//
//        pRes->pLeft = copy2(thiz->pLeft, leafList);
//        pRes->pRight = copy2(thiz->pRight, leafList);
//
//        if (pRes->pLeft) pRes->pLeft->Parent = pRes;
//        if (pRes->pRight) pRes->pRight->Parent = pRes;
//
//        return pRes;
//    }
//
//    static Relation CompressCSGTreeWithInside(CSGTree* tree, unsigned Id)
//    {
//        auto leaf = tree->Leaves[Id];
//
//        CSGTree::Node *curPtr = leaf, *parent = leaf->Parent;
//        while (parent && parent->Type == TYPE_UNION)
//        {
//            curPtr = parent;
//            parent = parent->Parent;
//        }
//
//        if (!parent)
//        {
//            delete tree->pRoot;
//            tree->Leaves.clear();
//            tree->pRoot = nullptr;
//            return REL_INSIDE;
//        }
//
//        assert(parent->Type == TYPE_INTERSECT);
//        int resCur = LeftOrRight(curPtr);
//        int resPar = LeftOrRight(parent);
//
//        if (resPar == 0) // root
//        {
//            if (resCur < 0)
//            {
//                tree->pRoot = parent->pRight;
//                parent->pRight = nullptr;
//                delete parent;
//            }
//            else
//            {
//                assert(resCur);
//                tree->pRoot = parent->pLeft;
//                parent->pLeft = nullptr;
//                delete parent;
//            }
//            tree->pRoot->Parent = nullptr;
//        }
//        else if (resPar < 0)
//        {
//            if (resCur < 0)
//            {
//                parent->pRight->Parent = parent->Parent;
//                parent->Parent->pLeft = parent->pRight;
//                parent->pRight = nullptr;
//                delete parent;
//            }
//            else
//            {
//                assert(resCur);
//                parent->pLeft->Parent = parent->Parent;
//                parent->Parent->pLeft = parent->pLeft;
//                parent->pLeft = nullptr;
//                delete parent;
//            }
//        }
//        else
//        {
//            if (resCur < 0)
//            {
//                parent->pRight->Parent = parent->Parent;
//                parent->Parent->pRight = parent->pRight;
//                parent->pRight = nullptr;
//                delete parent;
//            }
//            else
//            {
//                assert(resCur);
//                parent->pLeft->Parent = parent->Parent;
//                parent->Parent->pRight = parent->pLeft;
//                parent->pLeft = nullptr;
//                delete parent;
//            }
//        }
//        GetLeafList(tree);
//        return REL_UNKNOWN;
//    }
//
//    static Relation CompressCSGTreeWithOutside(CSGTree* tree, unsigned Id)
//    {
//        auto leaf = tree->Leaves[Id];
//        CSGTree::Node *curPtr = leaf, *parent = leaf->Parent;
//        while (parent && parent->Type == TYPE_INTERSECT)
//        {
//            curPtr = parent;
//            parent = parent->Parent;
//        }
//
//        if (!parent)
//        {
//            delete tree->pRoot;
//            tree->pRoot = nullptr;
//            tree->Leaves.clear();
//            return REL_OUTSIDE;
//        }
//
//        assert(parent->Type == TYPE_UNION);
//        int resCur = LeftOrRight(curPtr);
//        int resPar = LeftOrRight(parent);
//
//        if (resPar == 0) // root
//        {
//            if (resCur < 0)
//            {
//                tree->pRoot = parent->pRight;
//                parent->pRight = nullptr;
//                delete parent;
//            }
//            else
//            {
//                assert(resCur);
//                tree->pRoot = parent->pLeft;
//                parent->pLeft = nullptr;
//                delete parent;
//            }
//            tree->pRoot->Parent = nullptr;
//        }
//        else if (resPar < 0)
//        {
//            if (resCur < 0)
//            {
//                parent->pRight->Parent = parent->Parent;
//                parent->Parent->pLeft = parent->pRight;
//                parent->pRight = nullptr;
//                delete parent;
//            }
//            else
//            {
//                assert(resCur);
//                parent->pLeft->Parent = parent->Parent;
//                parent->Parent->pLeft = parent->pLeft;
//                parent->pLeft = nullptr;
//                delete parent;
//            }
//        }
//        else
//        {
//            if (resCur < 0)
//            {
//                parent->pRight->Parent = parent->Parent;
//                parent->Parent->pRight = parent->pRight;
//                parent->pRight = nullptr;
//                delete parent;
//            }
//            else
//            {
//                assert(resCur);
//                parent->pLeft->Parent = parent->Parent;
//                parent->Parent->pRight = parent->pLeft;
//                parent->pLeft = nullptr;
//                delete parent;
//            }
//        }
//        GetLeafList(tree);
//        return REL_UNKNOWN;
//    }
//
//    static Relation CompressCSGTreeWithSame(CSGTree* tree, unsigned Id)
//    {
//        auto leaf = tree->Leaves[Id];
//        CSGTree::Node *curPtr = leaf, *parent = leaf->Parent, *neib;
//
//        while (true)
//        {
//            if (!parent)
//            {
//                delete tree->pRoot;
//                tree->pRoot = nullptr;
//                tree->Leaves.clear();
//                return REL_SAME;
//            }
//            if (LeftOrRight(curPtr) < 0.0)
//                neib = parent->pRight;
//            else neib = parent->pLeft;
//
//            if (neib->relation == REL_SAME)
//            {
//                curPtr = parent;
//                parent = parent->Parent;
//            }
//            else if (neib->relation == REL_OPPOSITE)
//            {
//                if (parent->Type == TYPE_UNION)
//                    return CompressCSGTreeWithInside(tree, Id);
//                else if (parent->Type == TYPE_INTERSECT)
//                    return CompressCSGTreeWithOutside(tree, Id);
//                else assert(0);
//            }
//            else break;
//        }
//        curPtr->relation = REL_SAME;
//        SAFE_RELEASE(curPtr->pLeft);
//        SAFE_RELEASE(curPtr->pRight);
//        GetLeafList(tree);
//        return REL_NOT_AVAILABLE;
//    }
//
//    static Relation CompressCSGTreeWithOppo(CSGTree* tree, unsigned Id)
//    {
//        auto leaf = tree->Leaves[Id];
//        CSGTree::Node *curPtr = leaf, *parent = leaf->Parent, *neib;
//
//        while (true)
//        {
//            if (!parent)
//            {
//                delete tree->pRoot;
//                tree->pRoot = nullptr;
//                tree->Leaves.clear();
//                return REL_OPPOSITE;
//            }
//            if (LeftOrRight(curPtr) < 0.0)
//                neib = parent->pRight;
//            else neib = parent->pLeft;
//
//            if (neib->relation == REL_OPPOSITE)
//            {
//                curPtr = parent;
//                parent = parent->Parent;
//            }
//            else if (neib->relation == REL_SAME)
//            {
//                if (parent->Type == TYPE_UNION)
//                    return CompressCSGTreeWithInside(tree, Id);
//                else if (parent->Type == TYPE_INTERSECT)
//                    return CompressCSGTreeWithOutside(tree, Id);
//                else assert(0);
//            }
//            else break;
//        }
//        curPtr->relation = REL_OPPOSITE;
//        SAFE_RELEASE(curPtr->pLeft);
//        SAFE_RELEASE(curPtr->pRight);
//        GetLeafList(tree);
//        return REL_NOT_AVAILABLE;
//    }
//
//    Relation CompressCSGTree(CSGTree* tree, unsigned Id, Relation rel)
//    {
//        // T∩N⇒N,  F∪N⇒N, F∩N⇒F, and T∪N⇒T
//        if (tree->Leaves.find(Id) == tree->Leaves.end())
//            return REL_UNKNOWN;
//
//        switch (rel)
//        {
//        case REL_INSIDE:
//            return CompressCSGTreeWithInside(tree, Id);
//            break;
//        case REL_OUTSIDE:
//            return CompressCSGTreeWithOutside(tree, Id);
//            break;
//        case REL_SAME:
//            return CompressCSGTreeWithSame(tree, Id);
//            break;
//        case REL_OPPOSITE:
//            return CompressCSGTreeWithOppo(tree, Id);
//            break;
//        default:
//            assert(0);
//            break;
//        }
//
//        return REL_UNKNOWN;
//    }
//
//    Relation CompressCSGNodeIteration(CSGTree::Node*& root)
//    {
//        if (IsLeaf(root)) return root->relation;
//
//        Relation rRight, rLeft;
//        rRight = CompressCSGNodeIteration(root->pRight);
//
//        if (root->Type == TYPE_UNION)
//        {
//            if (rRight == REL_INSIDE)
//            {
//                delete root;
//                root = nullptr;
//                return REL_INSIDE;
//            }
//            else if (rRight == REL_OUTSIDE)
//            {
//                delete root->pRight;
//                auto tmp = root;
//                root = root->pLeft;
//                root->Parent = tmp->Parent;
//
//                tmp->pLeft = nullptr;
//                tmp->pRight = nullptr;
//                delete tmp;
//
//                return CompressCSGNodeIteration(root);
//            }
//        }
//#ifdef _DEBUG
//        else if (root->Type == TYPE_INTERSECT)
//#else
//        else
//#endif
//        {
//            if (rRight == REL_OUTSIDE)
//            {
//                delete root;
//                root = nullptr;
//                return REL_OUTSIDE;
//            }
//            else if (rRight == REL_INSIDE)
//            {
//                delete root->pRight;
//                auto tmp = root;
//                root = root->pLeft;
//                root->Parent = tmp->Parent;
//
//                tmp->pLeft = nullptr;
//                tmp->pRight = nullptr;
//                delete tmp;
//
//                return CompressCSGNodeIteration(root);
//            }
//        }
//#ifdef _DEBUG
//        else assert(0);
//#endif
//
//        rLeft = CompressCSGNodeIteration(root->pLeft);
//        if (root->Type == TYPE_UNION)
//        {
//            if (rLeft == REL_INSIDE)
//            {
//                delete root;
//                root = nullptr;
//                return REL_INSIDE;
//            }
//            else if (rLeft == REL_OUTSIDE)
//            {
//                delete root->pLeft;
//                auto tmp = root;
//                root = root->pRight;
//                root->Parent = tmp->Parent;
//
//                tmp->pLeft = nullptr;
//                tmp->pRight = nullptr;
//                delete tmp;
//                return rRight;
//            }
//        }
//#ifdef _DEBUG
//        else if (root->Type == TYPE_INTERSECT)
//#else
//        else
//#endif
//        {
//            if (rLeft == REL_OUTSIDE)
//            {
//                delete root;
//                root = nullptr;
//                return REL_OUTSIDE;
//            }
//            else if (rLeft == REL_INSIDE)
//            {
//                delete root->pLeft;
//                auto tmp = root;
//                root = root->pRight;
//                root->Parent = tmp->Parent;
//
//                tmp->pLeft = nullptr;
//                tmp->pRight = nullptr;
//                delete tmp;
//                return rRight;
//            }
//        }
//#ifdef _DEBUG
//        else assert(0);
//#endif
//
//        // TO DO: if no-one is outside or inside
//        // Suppose we do not have on and oppo relation
//        return REL_NOT_AVAILABLE;
//    }
//
//
//    Relation ParsingCSGTree(MPMesh* pMesh, Relation* tab, unsigned nMesh, CSGTree::Node* curTree, CSGTree::Node** leaves, TestTree& output)
//    {
//        for (unsigned i = 0; i < nMesh; i++)
//        {
//            if (leaves[i])
//                leaves[i]->relation = tab[i];
//        }
//
//        CSGTree::Node *seed = leaves[pMesh->ID], *comp;
//        int checkRel;
//        bool pass = true;
//        bool simple = true;
//        while (seed->Parent)
//        {
//            // T∩N⇒N,  F∪N⇒N, F∩N⇒F, and T∪N⇒T
//            if (seed->Parent->Type == TYPE_UNION)
//                checkRel = REL_OUTSIDE;
//            else
//#ifdef _DEBUG
//                if (seed->Parent->Type == TYPE_INTERSECT)
//#endif
//                    checkRel = REL_INSIDE;
//#ifdef _DEBUG
//                else assert(0);
//#endif
//                if (LeftOrRight(seed) < 0)
//                {
//                    comp = seed->Parent->pRight;
//                    checkRel ^= REL_SAME;
//                    seed->Parent->pRight = nullptr;
//                }
//                else
//                {
//                    comp = seed->Parent->pLeft;
//                    seed->Parent->pLeft = nullptr;
//                }
//                comp->Parent = nullptr;
//                Relation resRel = CompressCSGNodeIteration(comp);
//                if (resRel == REL_NOT_AVAILABLE)
//                {
//                    output.emplace_back();
//                    output.back().targetRelation = checkRel;
//                    output.back().testTree = comp;
//                    simple = false;
//                }
//                else
//                {
//                    SAFE_RELEASE(comp);
//                    if (!(checkRel & resRel))
//                    {
//                        pass = false;
//                        break;
//                    }
//                }
//
//                seed = seed->Parent;
//        }
//
//        SAFE_RELEASE(curTree);
//        if (pass)
//        {
//            if (!simple) return REL_NOT_AVAILABLE;
//            else return REL_SAME;
//        }
//        output.clear();
//        return REL_INSIDE; //也有可能是OutSide
//    }
//
//    CSGTree::Node* GetFirstNode(CSGTree::Node* root)
//    {
//        assert(root);
//        if (IsLeaf(root)) return root;
//        else	 return GetFirstNode(root->pLeft);
//    }
//
//    // In ∩ X = X;		Out ∩ X = Out;
//    // In ∪ X = In;		Out ∪ X = X;
//
//    CSGTree::Node* ParsingInside(CSGTree::Node*& curNode, Relation &output);
//    CSGTree::Node* ParsingOutside(CSGTree::Node*& curNode, Relation &output);
//    CSGTree::Node* ParsingSame(CSGTree::Node*& curNode, Relation &output);
//    CSGTree::Node* ParsingOppo(CSGTree::Node*& curNode, Relation &output);
//
//    CSGTree::Node* ParsingInside(CSGTree::Node*& curNode, Relation &output)
//    {
//        CSGTree::Node *neib, *parent = curNode->Parent;
//        double lor = LeftOrRight(curNode);
//        if (lor == 0)
//        {
//            output = REL_INSIDE;
//            return nullptr;
//        }
//
//        if (parent->Type == TYPE_INTERSECT)
//        {
//            if (lor < 0) // left node
//            {
//                neib = parent->pRight;
//                return GetFirstNode(neib);
//            }
//            else
//#ifdef _DEBUG
//                if (lor > 0) // right node, check the left
//#endif
//                {
//                    neib = parent->pLeft;
//                    if (neib->mark == mark) // same
//                        return ParsingSame(parent, output);
//                    else if (neib->mark == mark + 1)
//                        return ParsingOppo(parent, output);
//                    else return ParsingInside(parent, output);
//                }
//#ifdef _DEBUG
//                else assert(0);
//#endif
//        }
//        else // Union
//            return ParsingInside(parent, output);
//    }
//
//    CSGTree::Node* ParsingOutside(CSGTree::Node*& curNode, Relation &output)
//    {
//        CSGTree::Node *neib, *parent = curNode->Parent;
//        double lor = LeftOrRight(curNode);
//        if (lor == 0)
//        {
//            output = REL_OUTSIDE;
//            return nullptr;
//        }
//
//        if (parent->Type == TYPE_UNION)
//        {
//            if (lor < 0) // left node
//            {
//                neib = parent->pRight;
//                return GetFirstNode(neib);
//            }
//            else
//#ifdef _DEBUG
//                if (lor > 0) // right node, check the left
//#endif
//                {
//                    neib = parent->pLeft;
//                    if (neib->mark == mark) // same
//                        return ParsingSame(parent, output);
//                    else if (neib->mark == mark + 1)
//                        return ParsingOppo(parent, output);
//                    else return ParsingOutside(parent, output);
//                }
//#ifdef _DEBUG
//                else assert(0);
//#endif
//        }
//        else // Union
//            return ParsingOutside(parent, output);
//    }
//
//    CSGTree::Node* ParsingSame(CSGTree::Node*& curNode, Relation &output)
//    {
//        CSGTree::Node *neib, *parent = curNode->Parent;
//        double lor = LeftOrRight(curNode);
//        if (lor == 0)
//        {
//            output = REL_SAME;
//            return nullptr;
//        }
//
//        if (lor < 0) // left
//        {
//            curNode->mark = mark;
//            return GetFirstNode(parent->pRight);
//        }
//        else
//        {
//            neib = parent->pLeft;
//            if (neib->mark == mark) // same
//                return ParsingSame(parent, output);
//            else if (neib->mark == mark + 1)
//            {
//                if (parent->Type == TYPE_INTERSECT)
//                    return ParsingOutside(parent, output);
//                else return ParsingInside(parent, output);
//            }
//            else return ParsingSame(parent, output);
//        }
//    }
//
//    CSGTree::Node* ParsingOppo(CSGTree::Node*& curNode, Relation &output)
//    {
//        CSGTree::Node *neib, *parent = curNode->Parent;
//        double lor = LeftOrRight(curNode);
//        if (lor == 0)
//        {
//            output = REL_OPPOSITE;
//            return nullptr;
//        }
//
//        if (lor < 0) // left
//        {
//            curNode->mark = mark + 1;
//            return GetFirstNode(parent->pRight);
//        }
//        else
//        {
//            neib = parent->pLeft;
//            if (neib->mark == mark + 1) // oppo
//                return ParsingOppo(parent, output);
//            else if (neib->mark == mark)
//            {
//                if (parent->Type == TYPE_INTERSECT)
//                    return ParsingOutside(parent, output);
//                else return ParsingInside(parent, output);
//            }
//            else return ParsingOppo(parent, output);
//        }
//    }
//
//    CSGTree::Node* GetNextNode(CSGTree::Node* curNode, Relation rel, Relation &output)
//    {
//        assert(curNode->Type == TYPE_LEAF);
//        CSGTree::Node *parent = curNode->Parent;
//
//        switch (rel)
//        {
//        case REL_INSIDE:
//            return ParsingInside(curNode, output);
//        case REL_OUTSIDE:
//            return ParsingOutside(curNode, output);
//        case REL_SAME:
//            return ParsingSame(curNode, output);
//        case REL_OPPOSITE:
//            return ParsingOppo(curNode, output);
//        default: break;
//        }
//        assert(0);
//    }

}  // namespace CSG

