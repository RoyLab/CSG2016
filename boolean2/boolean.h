#pragma once
#include <string>
#include <vector>
#include "macros.h"
#include "global.h"

extern "C"
{
	void XRWY_DLL solveCSG(const std::string& expr, std::vector<Boolean::RegularMesh*>& meshes);
	void XRWY_DLL test(std::vector<std::string>& names, std::string& expr, std::string& output);
}

