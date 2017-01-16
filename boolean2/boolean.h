#pragma once
#include <string>
#include <vector>
#include <macros.h>
#include "global.h"

struct CsgInputs
{
    std::vector<std::string> names;
    std::string expr;
    std::string output;
    std::string log;

    bool need_log = false;
    bool need_pwn_test = false;
};

extern "C"
{

	XRWY_DLL Boolean::RegularMesh* solveCSG(const std::string& expr, std::vector<Boolean::RegularMesh*>& meshes);
	XRWY_DLL void test(const CsgInputs& inputs);

    XRWY_DLL void test1();
    XRWY_DLL void test2();
    XRWY_DLL void test3();
    XRWY_DLL void test4();
}

