#pragma once
#include <macros.h>
#include <vector>
#include <string>



extern "C" XRWY_DLL void cgaleval(std::vector<std::string>& names, const std::string& expr, const std::string& output, const std::string& log);
extern "C" XRWY_DLL void corkeval(std::vector<std::string>& names, const std::string& expr, const std::string& output, const std::string& log);
extern "C" XRWY_DLL void cgaleval_test();
