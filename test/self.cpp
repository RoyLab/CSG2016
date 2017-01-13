#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>

#include <cmdline.h>

#include "boolean.h"
#include "cgaleval.h"


void TestByMethod(std::vector<std::string>& names, std::string& expr, 
    const std::string& output, std::string method, const std::string& tmp_log);

int self_main(int argc, char* argv[])
{
    cmdline::parser cmd_parser;
    cmd_parser.add<std::string>("name", 'n', "self intersection", true);
    cmd_parser.add<std::string>("method", 'm', "method used for evaluation", false, "my", cmdline::oneof<std::string>("cgal", "my", "cork"));
    cmd_parser.add<std::string>("output", 'o', "output file name", true);
    cmd_parser.add<std::string>("tmplop", 'l', "log to record", false, "");

    cmd_parser.parse_check(argc, argv);
    std::string filename = cmd_parser.get<std::string>("name");
    std::string method = cmd_parser.get<std::string>("method");
    std::string output = cmd_parser.get<std::string>("output");
    std::string temp_log = cmd_parser.get<std::string>("tmplop");

    std::string expr;
    std::vector<std::string> names;

    expr = "0+1";
    names.clear();
    names.push_back(filename);
    names.push_back(filename);

    try
    {
        TestByMethod(names, expr, output, method, temp_log);
    }
    catch (...)
    {
        std::cerr << "IO error\n";
        system("pause");
        return 1;
    }

    system("pause");
    return 0;
}