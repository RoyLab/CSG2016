#pragma once
#include <macros.h>
#include <vector>
#include <string>



extern "C" XRWY_DLL void cgaleval(std::vector<std::string>& names, std::string& expr, const std::string& output);
extern "C" XRWY_DLL void cgaleval_test();
