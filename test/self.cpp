#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>

#include <cmdline.h>

#include "boolean.h"
#include "cgaleval.h"


void TestByMethod(std::vector<std::string>& names, std::string& expr, const std::string& output, std::string method);

int main2(int argc, char* argv[])
{
    cmdline::parser cmd_parser;
    cmd_parser.add<std::string>("name", 'n', "self intersection", true);
    cmd_parser.add<std::string>("method", 'm', "method used for evaluation", false, "my", cmdline::oneof<std::string>("cgal", "my", "cork"));
    cmd_parser.add<std::string>("output", 'o', "output file name", true);

    cmd_parser.parse_check(argc, argv);
    std::string filename = cmd_parser.get<std::string>("name");
    std::string method = cmd_parser.get<std::string>("method");
    std::string output = cmd_parser.get<std::string>("output");

    std::string expr;
    std::vector<std::string> names;

    expr = "0+1";
    names.clear();
    names.push_back(filename);
    names.push_back(filename);
    TestByMethod(names, expr, output, method);

    return 0;
}