#define _USE_MATH_DEFINES
#include "bool.h"
#include <iostream>
using namespace CSG;

int main()
{
    std::vector<std::string> names;
    names.push_back("../../models/box1.off");
    names.push_back("../../models/box4.off");

    std::string expr("0*1");
    test(names, expr);
    std::cout << "This is the end!" << std::endl;
    system("pause");
    return 0;
}