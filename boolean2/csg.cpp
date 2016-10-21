#include "precompile.h"
#include "csg.h"
#include <stack>

namespace Boolean
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

} 

