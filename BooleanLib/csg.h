#pragma once
//#include "COctree.h"
#include "macroutil.h"
#include <stack>
#include <map>
#include <vector>

namespace CSG
{

    enum Relation
    {
        REL_UNKNOWN =   0x0000,
        REL_INSIDE =    0x0001,
        REL_OUTSIDE =   0x0002,
        REL_SAME =      0x0004,
        REL_OPPOSITE =  0x0008,

        REL_NOT_AVAILABLE = -1 // (0xffffffff)
    };

    enum CSGNodeOp
    {
        TYPE_UNKNOWN = 0,
        TYPE_UNION, TYPE_INTERSECT, TYPE_DIFF, TYPE_LEAF
    };



    static inline bool IsOperator(char c)
    {
        if ((c == '(') || (c == ')') || (c == '+')
            || (c == '*') || (c == '-'))
            return true;
        return false;
    }


    std::string InfixToPostfix(const std::string& infix);



    template <class Mesh>
    class CSGTree
    {
        template <class Mesh>
        struct NodeBase
        {
            CSGNodeOp type = TYPE_UNKNOWN;
            Relation relation = REL_UNKNOWN;
            NodeBase *pLeft = nullptr, *pRight = nullptr, *Parent = nullptr;
            //int mark = -1;

            Mesh* pMesh = nullptr;
            bool bInverse = false;

            virtual ~NodeBase(){}
        };

        typedef NodeBase<Mesh> TreeNode;

    public:

        struct LeafNode:
            public TreeNode
        {
            LeafNode(Mesh* mesh) : TreeNode()
            {
                TreeNode::type = TYPE_LEAF;
                TreeNode::pMesh = mesh;
            }
        };

        struct OpNode:
            public TreeNode
        {
            OpNode(CSGNodeOp op, TreeNode* left, TreeNode* right) : TreeNode()
            {
                TreeNode::type = op;
                TreeNode::pLeft = left;
                TreeNode::pRight = right;
            }
        };

        CSGTree() : pRoot(nullptr) {}
        CSGTree(const CSGTree& other) { pRoot = copyNodeWithoutMesh(other.pRoot); }
        ~CSGTree(){ reset(); }

        void reset()
        {
            for (TreeNode* node : nodeList)
                delete node;
            nodeList.clear();

            pRoot = nullptr;
        }

        TreeNode* copyNodeWithoutMesh(const TreeNode& root)
        {
            TreeNode *res = reg(new TreeNode(root));
            TreeNode* left, right;

            if (root.pLeft && root.pLeft)
            {
                left = copyNode(root.pLeft);
                right = copyNode(root.pRight);

                res->pLeft = left;
                res->pRight = right;
                left->pParent = res;
                right->pParent = res;
            }

            return res;
        }

        void createCSGTreeFromExpr(const std::string& expr, Mesh** meshes, size_t nMesh)
        {
            std::string postfix = InfixToPostfix(expr);
            std::stack<TreeNode*> operandstack;
            int i = 0;

            while (i < postfix.size())
            {
                int j = 0;
                char number[16];
                while (i< postfix.size())
                {
                    if (postfix[i] == ' ' || IsOperator(postfix[i]))
                    {
                        if (j)
                        {
                            number[j] = 0;
                            int idx = atoi(number);
                            assert(idx < nMesh);
                            operandstack.push(reg(new LeafNode(meshes[idx])));
                            j = 0;
                        }
                        if (IsOperator(postfix[i]))
                            break;
                        i++;
                    }
                    else number[j++] = postfix[i++];
                }
                if (i < postfix.size())
                {
                    TreeNode* right = operandstack.top();
                    operandstack.pop();
                    TreeNode* left = operandstack.top();
                    operandstack.pop();
                    switch (postfix[i])
                    {
                    case '+':
                        operandstack.push(reg(new OpNode(TYPE_UNION, left, right)));
                        break;
                    case '*':
                        operandstack.push(reg(new OpNode(TYPE_INTERSECT, left, right)));
                        break;
                    case '-':
                        operandstack.push(reg(new OpNode(TYPE_DIFF, left, right)));
                        break;
                    default:
                        assert(0);
                    }
                    i++;
                }


            }
            pRoot = operandstack.top();
            operandstack.pop();
        }

        void makePositiveAndLeftHeavy()
        {
            size_t maxLvl = 0;
            pRoot = makePositiveAndLeftHeavy(pRoot, false, 0u, maxLvl);
        }

        TreeNode* makePositiveAndLeftHeavy(TreeNode* root, bool inverse, size_t lvl, size_t &maxLvl)
        {
            if (root->type == TYPE_LEAF)
            {
                root->bInverse = inverse;
                maxLvl = (maxLvl < lvl) ? lvl : maxLvl;
            }
            else
            {
                size_t Ldepth(0), Rdepth(0);

                if (root->type == TYPE_DIFF)
                {
                    root->type = TYPE_INTERSECT;
                    makePositiveAndLeftHeavy(root->pLeft, inverse, lvl + 1, Ldepth);
                    makePositiveAndLeftHeavy(root->pRight, !inverse, lvl + 1, Rdepth);
                }
                else
                {
                    makePositiveAndLeftHeavy(root->pLeft, inverse, lvl + 1, Ldepth);
                    makePositiveAndLeftHeavy(root->pRight, inverse, lvl + 1, Rdepth);
                }

                maxLvl = (Ldepth > Rdepth) ? Ldepth : Rdepth;

                if (Ldepth < Rdepth)
                {
                    auto tmpNode = root->pLeft;
                    root->pLeft = root->pRight;
                    root->pRight = tmpNode;
                }

                if (inverse)
                    root->type = (root->type == TYPE_INTERSECT) ? TYPE_UNION : TYPE_INTERSECT;
            }
            
            return root;
        }

    private:
        TreeNode* reg(TreeNode* node)
        {
            nodeList.push_back(node);
            return node;
        }

        TreeNode* pRoot;
        std::vector<TreeNode*> nodeList;
    };




    //struct Branch
    //{
    //    int targetRelation;
    //    Node* testTree;

    //    ~Branch() { SAFE_DELETE(testTree); }
    //};

    //typedef std::list<Branch> TestTree;
    //void GetLeafList(CSGTreeNode* root, std::vector<int>& list);
    //CSGTree* ConvertCSGTree(GS::CSGExprNode* root, MPMesh*** arrMesh, int *nMes); // convert nodes.
    //CSGTree* ConvertToPositiveTree(const CSGTree* tree);
    //Relation CompressCSGTree(CSGTree* tree, unsigned Id, Relation rel);
    //Relation ParsingCSGTree(MPMesh* pMesh, Relation* tab, unsigned nMesh, CSGTreeNode* curTree, CSGTreeNode** leaves, TestTree& output);
    //CSGTreeNode* GetNextNode(CSGTreeNode* curNode, Relation rel, Relation &output);
    //CSGTreeNode* GetFirstNode(CSGTreeNode* root);
    //inline bool IsLeaf(CSGTreeNode* node) { return !(node->pLeft && node->pRight); }
    //CSGTree* copy(const CSGTree* thiz);
    //CSGTreeNode* copy2(const CSGTreeNode* thiz, CSGTreeNode** leafList);
    //Relation CompressCSGNodeIteration(CSGTreeNode*& root);

    ///**
    //if it is a left child, return negative
    //if it is a right child, return positive
    //if it is a root, return 0
    //*/
    //inline int LeftOrRight(CSGTreeNode* node)
    //{
    //    assert(node);

    //    if (!node->Parent) return 0;
    //    if (node->Parent->pLeft == node) return -1;
    //    if (node->Parent->pRight == node) return 1;

    //    assert(0);
    //    return 0;
    //}

} // namespace CSG

