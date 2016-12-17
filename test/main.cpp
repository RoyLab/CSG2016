#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>

#include <cmdline.h>

#include "boolean.h"
#include "cgaleval.h"

std::stringstream getValidStream(std::ifstream& file)
{
    std::string line;
    while (file)
    {
        std::getline(file, line);
        if (line.size() && line[0] != ';')
            return (std::stringstream(line));
    }
    throw std::exception();
}

void TestByMethod(std::vector<std::string>& names, std::string& expr, const std::string& output, std::string method)
{
    if (method == "my")
        test(names, expr, output);
    else if (method == "cgal")
        cgaleval(names, expr, output);
    else if (method == "cork")
        corkeval(names, expr, output);
}

int main(int argc, char* argv[])
{
    cmdline::parser cmd_parser;
    cmd_parser.add<std::string>("script", 's', "the script to run", false, "./mycsg.ini");
    cmd_parser.add<std::string>("method", 'm', "method used for evaluation", false, "my", cmdline::oneof<std::string>("cgal", "my", "cork"));

    cmd_parser.parse_check(argc, argv);
    std::string config_filename = cmd_parser.get<std::string>("script");
    std::string method = cmd_parser.get<std::string>("method");

    std::string expr;
    std::vector<std::string> names;

    if (false && argc == 1)
    {
        expr = "0+1";
        names.clear();
        names.push_back("D:/Codes/Boolean2016/exps/data/cmp_meshworks/buddha.off");
        names.push_back("D:/Codes/Boolean2016/exps/data/cmp_meshworks/lion.off");
        TestByMethod(names, expr, "D:/result.off", method);
        system("pause");
        return 0;
    }

    std::ifstream configFile(config_filename);
    if (!configFile.is_open())
    {
        std::cerr << "IO error: " << config_filename << std::endl;
        return -1;
    }
    int nMesh;
    std::string name, output, line;
    std::stringstream buffer;
    
    //try
    //{
        while (1)
        {
            nMesh = 0;
            names.clear();
            buffer = getValidStream(configFile);
            buffer >> nMesh;
            if (nMesh == 0 || configFile.eof()) break;

            for (int j = 0; j < nMesh; j++)
            {
                buffer = getValidStream(configFile);
                buffer >> name;
                if (name[0] == '@')
                {
                    --j;
                    name = name.substr(1, name.length());
                    buffer = getValidStream(configFile);
                    int s, e;
                    char nameBuffer[512];
                    buffer >> s >> e;
                    for (int k = s; k <= e; k++, j++)
                    {
                        sprintf(nameBuffer, name.c_str(), k);
                        names.push_back(nameBuffer);
                    }
                }
                else names.push_back(name);
            }
            // expr
            buffer = getValidStream(configFile);
            char first = buffer.peek();
            if (first == '@')
            {
                buffer.get();
                expr = "0";
                int num, count = 1;
                std::string opstr;
                char opch;
                while (1)
                {
                    buffer >> num >> opstr;
                    assert(opstr.size() == 1);
                    opch = opstr[0];
                    for (int j = 0; j < num; j++, count++)
                    {
                        expr += opch;
                        expr += std::to_string(count);
                    }
                    if (buffer.eof()) break;
                }
            }
            else buffer >> expr;

            // output
            buffer = getValidStream(configFile);
            buffer >> output;

            TestByMethod(names, expr, output, method);
        }
    //}
    //catch (...)
    //{
    //    std::cerr << "non-standard ends or error occurred.\n";
    //}

    configFile.close();

    if (argc == 1)
        system("pause");
    return 0;
}