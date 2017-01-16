#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>

#include <cmdline.h>

#include "boolean.h"
#include "cgaleval.h"

void TestSelfInsct(const std::string& name, const std::string& output,
    const std::string& method, const std::string& temp)
{
    std::vector<std::string> names(2, name);
    std::string expr = "0+1";
    if (method == "my")
    {
        CsgInputs inputs;
        inputs.names = names;
        inputs.expr = expr;
        inputs.output = output;
        inputs.log = temp;

        inputs.need_log = true;
        inputs.need_pwn_test = true;
        test(inputs);
    }
    else if (method == "cgal")
        cgaleval(names, expr, output, temp);
    else if (method == "cork")
        corkeval(names, expr, output, temp);
}


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

    try
    {
        TestSelfInsct(filename, output, method, temp_log);
    }
    catch (...)
    {
        std::cerr << "IO error\n";
        return 1;
    }

    return 0;
}