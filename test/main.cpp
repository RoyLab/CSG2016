#define _USE_MATH_DEFINES
#include "bool.h"
#include <iostream>
using namespace CSG;

int main()
{

	std::vector<int> a;
	int b[20] = { 2,3, 4,5, 6, 7 };

	a.insert(a.end(), b, b + 20);

	for (int i : a)
		std::cout << i << std::endl;
	system("pause");
	return 0;
    //std::vector<std::string> names;
    //names.push_back("../../models/bunny.off");
    //names.push_back("../../models/dragon.off");

    //std::string expr("0+1");
    //test(names, expr);
    //std::cout << "This is the end!" << std::endl;
    //system("pause");
    //return 0;
}