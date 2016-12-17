#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <cassert>
#include <cctype>
using std::string;
using std::vector;

#include <cmdline.h>

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

int main(int argc, char* argv[])
{
    cmdline::parser cmd_parser;
    cmd_parser.add<std::string>("script", 's', "the script to run", false, "../test/mycsg.ini");
    cmd_parser.add<std::string>("output", 'o', "the outputfile", false, "D:/test.bat");

    cmd_parser.parse_check(argc, argv);
    std::string config_filename = cmd_parser.get<std::string>("script");
    std::string output_filename = cmd_parser.get<std::string>("output");

    string expr;
    vector<string> names;
    std::ifstream configFile(config_filename);
    std::ofstream outputfile(output_filename, 'w');
    //auto &outputfile = std::cout;
    if (!configFile.is_open())
    {
        std::cerr << "IO error: " << config_filename << std::endl;
        return -1;
    }
    int nMesh;
    std::string name, output, line;
    std::stringstream buffer;

    try
    {
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

            for (string &name : names)
            {
                if (name[1] == ':')
                    name = name.substr(2);
                name.replace(name.rfind("off"), 3, "obj");
            }
            if (output[1] == ':')
                output = output.substr(2);
            output.replace(output.rfind("off"), 3, "obj");

            outputfile << "carve -O \"";
            int mesh_number = 0;
            bool number_flag = false;
            for (char ch : expr)
            {
                if (std::isdigit(ch))
                {
                    if (!number_flag)
                    {
                        number_flag = true;
                    }
                }
                else
                {
                    if (number_flag)
                    {
                        number_flag = false;
                        string name = names[mesh_number++];
                        outputfile << name;
                    }

                    switch (ch)
                    {
                    case '+':
                        outputfile << " | "; break;
                    case '*':
                        outputfile << " & "; break;
                    default:
                        outputfile << ' ' << ch << ' ';
                        break;
                    }
                }

            }
            if (mesh_number < names.size())
            {
                assert(mesh_number == names.size() - 1);
                outputfile << names.back();
            }
            outputfile << "\" > " << output;
            outputfile << std::endl;
        }
    }
    catch (...)
    {
        std::cerr << "non-standard ends or error occurred.\n";
    }

    configFile.close();
    outputfile.close();

    //system("pause");
    return 0;
}
