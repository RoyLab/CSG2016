#define _USE_MATH_DEFINES
#include <iostream>
#include "boolean.h"

int main()
{
    std::vector<std::string> names;
    names.push_back("../../models/ball1.off");
    names.push_back("../../models/ball2.off");

    std::string expr("0+1");
	test2(names, expr, expr);
    //std::cout << "This is the end!" << std::endl;
    system("pause");
    return 0;
}