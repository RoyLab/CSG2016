#define _USE_MATH_DEFINES
#include <iostream>
#include "boolean.h"

int main()
{
    //names.push_back("../../models/box1.off");
    //names.push_back("../../models/box2.off");

    //names.push_back("../../models/ball1.off");
    //names.push_back("../../models/ball2.off");
    std::string expr;
    std::vector<std::string> names;
    if (true)
    {
        expr = "0+1";
        names.clear();
        names.push_back("../../models/box1.off");
        names.push_back("../../models/box2.off");
        test(names, expr, "D:/result1.off");

        names.clear();
        names.push_back("../../models/box1.off");
        names.push_back("D:/offlib/box3.off");
        test(names, expr, "D:/result2.off");

        names.clear();
        names.push_back("../../models/ball1.off");
        names.push_back("../../models/ball2.off");
        test(names, expr, "D:/result3.off");

        names.clear();
        names.push_back("../../models/ball1.off");
        names.push_back("../../models/ball1.off");
        test(names, expr, "D:/result4.off");

        names.clear();
        names.push_back("../../models/box-lean-1.off");
        names.push_back("../../models/box-lean-2.off");
        test(names, expr, "D:/result5.off");

#ifndef _DEBUG
        names.clear();
        names.push_back("../../models/bunny.off");
        names.push_back("../../models/dragon.off");
        test(names, expr, "D:/result6.off");
#endif // !_DEBUG

        expr = "0+1+2";
        names.clear();
        names.push_back("../../models/box1.off");
        names.push_back("../../models/box2.off");
        names.push_back("../../models/boxm1.off");
        test(names, expr, "D:/result7.off");

        expr = "0+1+2";
        names.clear();
        names.push_back("../../models/box1.off");
        names.push_back("../../models/box2.off");
        names.push_back("../../models/box3.off");
        test(names, expr, "D:/result8.off");

        expr = "0*1";
        names.clear();
        names.push_back("../../models/box1.off");
        names.push_back("../../models/tetrahedron.off");
        test(names, expr, "D:/result9.off");

    }

    expr = "0";
    char buffer[200];
    for (int i = 1; i < 25; i++)
    {
        expr += "+";
        expr += itoa(i, buffer, 10);
    }

    for (int i = 25; i < 50; i++)
    {
        expr += "-";
        expr += itoa(i, buffer, 10);
    }

    names.clear();
    for (int i = 0; i < 50; i++)
    {
        sprintf(buffer, "D:/Codes/Boolean2016/exps/data/ref_timing/t1_%d.off", i);
        names.push_back(buffer);
    }
    test(names, expr, "D:/result.off");

    system("pause");
    return 0;
}