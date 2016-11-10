#define _USE_MATH_DEFINES
#include <iostream>
#include "boolean.h"

int main()
{
    //names.push_back("../../models/box1.off");
    //names.push_back("../../models/box2.off");

    //names.push_back("../../models/ball1.off");
    //names.push_back("../../models/ball2.off");

    std::string expr("0-1");
    std::vector<std::string> names;
    if (false)
    {
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

    }
    names.clear();
    names.push_back("../../models/bunny.off");
    names.push_back("../../models/dragon.off");
    test(names, expr, "D:/result.off");

    system("pause");
    return 0;
}