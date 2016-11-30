#pragma once
#include <string>
#include <vector>
#include <stack>
#include <cassert>

namespace boolean
{
    enum CsgOperator
    {
        CSG_OP_UNKNOWN = 0,
        CSG_OP_UNION, 
        CSG_OP_INTERSECT,
        CSG_OP_DIFF
    };

    template <class MeshType>
    struct NodeBase
    {
        CsgOperator type;
        MeshType    *meshRef;
        NodeBase    *left, *right;
        bool        is_inverse;
    };
    
    template <class MeshType>
    class CsgBinaryEval
    {
    public:
        virtual MeshType* operator()(MeshType* a, MeshType* b, CsgOperator op) const = 0;
    };

    template <class MeshType>
    class CsgTree
    {
    private:
        typedef NodeBase<MeshType> TreeNode;
    private:
        TreeNode* root_;
        std::vector<TreeNode*> nodelist_;

    public:
        CsgTree():root_(nullptr) {}
        ~CsgTree();

        void createCSGTreeFromExpr(const std::string& expr, std::vector<MeshType*>& meshes);
        MeshType* eval(CsgBinaryEval<MeshType>* method) { return eval(root_, method); }
    private:
        MeshType* eval(TreeNode* node, CsgBinaryEval<MeshType>* method);
        TreeNode* reg(TreeNode* node);
    };

    /***************************************************************************************************/
    /***************************************************************************************************/
    /***********************************Implementation**************************************************/
    /***************************************************************************************************/
    /***************************************************************************************************/
    namespace
    {
        bool IsOperator(char c)
        {
            if ((c == '(') || (c == ')') || (c == '+')
                || (c == '*') || (c == '-'))
                return true;
            return false;
        }

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
    }

    template <class MeshType>
    inline CsgTree<MeshType>::~CsgTree()
    {
        for (TreeNode* node : nodelist_)
            delete node;
    }

    template <class MeshType>
    inline void CsgTree<MeshType>::createCSGTreeFromExpr(const std::string& expr, std::vector<MeshType*>& meshes)
    {
        std::string postfix = InfixToPostfix(expr);
        std::stack<TreeNode*> operandstack;
        int i = 0;
        size_t nMesh = meshes.size();
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
                        operandstack.push(reg(new TreeNode{ CSG_OP_UNKNOWN, meshes[idx] }));
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
                    operandstack.push(reg(new TreeNode{ CSG_OP_UNION, 0, left, right }));
                    break;
                case '*':
                    operandstack.push(reg(new TreeNode{ CSG_OP_INTERSECT, 0, left, right }));
                    break;
                case '-':
                    operandstack.push(reg(new TreeNode{ CSG_OP_DIFF, 0, left, right }));
                    break;
                default:
                    assert(0);
                }
                i++;
            }


        }
        root_ = operandstack.top();
        operandstack.pop();
    }

    template<class MeshType>
    inline MeshType * CsgTree<MeshType>::eval(TreeNode * node, CsgBinaryEval<MeshType>* method)
    {
        if (node->type == CSG_OP_UNKNOWN)
        {
            return node->meshRef;
        }
        else
        {
            MeshType* left_eval = eval(node->left, method);
            MeshType* right_eval = eval(node->right, method);
            MeshType* result = (*method)(left_eval, right_eval, node->type);

            if (node->left->type != CSG_OP_UNKNOWN)
            {
                delete left_eval;
            }

            if (node->right->type != CSG_OP_UNKNOWN)
            {
                delete right_eval;
            }
            return result;
        }
    }

    template <class MeshType>
    inline typename CsgTree<MeshType>::TreeNode* CsgTree<MeshType>::reg(TreeNode* node)
    {
        nodelist_.push_back(node);
        return node;
    }
}
