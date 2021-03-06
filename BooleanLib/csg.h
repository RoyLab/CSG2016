#pragma once
#include <stack>
#include <map>
#include <vector>
#include <cstdint>
#include <unordered_map>

#include "macroutil.h"
#include "BinaryTree.h"

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

    typedef int8_t Indicator;

    class IndicatorVector
    {
    public:
        virtual Indicator& operator[](size_t meshId) = 0;
        virtual const Indicator& operator[](size_t meshId) const = 0;
    };


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
                left->Parent = this;
                right->Parent = this;
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
            meshNodeList.clear();

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
            meshNodeList.resize(nMesh);

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
                root->pMesh->bInverse = inverse;
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

        Relation eval(int meshId, IndicatorVector& inds)
        {
            assert(oldTree);

        }

        TreeNode* getRoot(){ return pRoot; }

    private:
        TreeNode* reg(TreeNode* node)
        {
            nodeList.push_back(node);
            if (node->type == TYPE_LEAF)
                meshNodeList[node->pMesh->Id] = node;
            return node;
        }

        Relation eval(Relation a, Relation b, CSGNodeOp op)
        {
            assert(a != REL_NOT_AVAILABLE);
            assert(b != REL_NOT_AVAILABLE);
            assert(op != TYPE_UNKNOWN);
            assert(op != TYPE_LEAF);
            assert(op != TYPE_DIFF);

            switch (op)
            {
            case TYPE_UNION:
                evalUnion(a, b);
                break;
            case TYPE_INTERSECT:
                evalIntersect(a, b);
                break;
            default:
                assert(0);
                break;
            }
        }

        Relation evalUnion(Relation a, Relation b)
        {
            assert(a != REL_UNKNOWN);
            assert(b != REL_UNKNOWN);

            if (a == REL_INSIDE || b == REL_INSIDE)
                return REL_INSIDE;

            if (a == REL_OUTSIDE)
                return b;

            if (b == REL_OUTSIDE)
                return a;

            if (a == b)
                return a;

            else return REL_INSIDE;
        }

        Relation evalIntersect(Relation a, Relation b)
        {
            assert(a != REL_UNKNOWN);
            assert(b != REL_UNKNOWN);

            if (a == REL_OUTSIDE || b == REL_OUTSIDE)
                return REL_OUTSIDE;

            if (a == REL_INSIDE)
                return b;

            if (b == REL_INSIDE)
                return a;

            if (a == b)
                return a;

            else return REL_OUTSIDE;
        }

        TreeNode* pRoot = nullptr;
        CSGTreeOld* oldTree = nullptr;
        std::vector<TreeNode*> nodeList, meshNodeList;
    };


} // namespace CSG

