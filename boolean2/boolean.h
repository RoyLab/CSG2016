#pragma once
#include <string>
#include <vector>
#include <macros.h>
#include "global.h"

extern "C"
{
	XRWY_DLL Boolean::RegularMesh* solveCSG(const std::string& expr, std::vector<Boolean::RegularMesh*>& meshes);
	XRWY_DLL void test(std::vector<std::string>& names, std::string& expr, const std::string& output);
}

