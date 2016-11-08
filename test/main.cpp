#define _USE_MATH_DEFINES
#include <iostream>
#include "boolean.h"

int main()
{
    std::vector<std::string> names;
    names.push_back("../../models/dragon.off");
    names.push_back("../../models/bunny.off");

    std::string expr("0-1");

	test(names, expr, "D:/result.off");
    system("pause");
    return 0;
}