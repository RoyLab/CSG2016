#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>
#include "boolean.h"

std::stringstream getValidStream(std::ifstream& file)
{
    std::string line;
    while (file)
    {
        std::getline(file, line);
        if (line.size() && line[0] != ';')
            return (std::stringstream(line));
    }
}

int main()
{
    std::string expr;
    std::vector<std::string> names;

    std::ifstream configFile(".\\mycsg.ini");
    //std::ifstream configFile("D:\\Codes\\Boolean2016\\exps\\mycsg.ini");
    if (!configFile.is_open())
    {
        if (true)
        {
            expr = "0+1";
            names.clear();
            names.push_back("../../models/box1.off");
            names.push_back("../../models/box2.off");
            test(names, expr, "D:/result1.off");

            names.clear();
            names.push_back("../../models/box1.off");
            names.push_back("D:/offlib/box3.off");
            test(names, expr, "D:/result2.off");

            names.clear();
            names.push_back("../../models/ball1.off");
            names.push_back("../../models/ball2.off");
            test(names, expr, "D:/result3.off");

            names.clear();
            names.push_back("../../models/ball1.off");
            names.push_back("../../models/ball1.off");
            test(names, expr, "D:/result4.off");

            names.clear();
            names.push_back("../../models/box-lean-1.off");
            names.push_back("../../models/box-lean-2.off");
            test(names, expr, "D:/result5.off");

#ifndef _DEBUG
            names.clear();
            names.push_back("../../models/bunny.off");
            names.push_back("../../models/dragon.off");
            test(names, expr, "D:/result6.off");
#endif // !_DEBUG

            expr = "0+1+2";
            names.clear();
            names.push_back("../../models/box1.off");
            names.push_back("../../models/box2.off");
            names.push_back("../../models/boxm1.off");
            test(names, expr, "D:/result7.off");

            expr = "0+1+2";
            names.clear();
            names.push_back("../../models/box1.off");
            names.push_back("../../models/box2.off");
            names.push_back("../../models/box3.off");
            test(names, expr, "D:/result8.off");

            expr = "0*1";
            names.clear();
            names.push_back("../../models/box1.off");
            names.push_back("../../models/tetrahedron.off");
            test(names, expr, "D:/result9.off");

            expr = "0";
            char buffer[200];
            for (int i = 1; i < 25; i++)
            {
                expr += "+";
                expr += itoa(i, buffer, 10);
            }

            for (int i = 25; i < 50; i++)
            {
                expr += "-";
                expr += itoa(i, buffer, 10);
            }

            names.clear();
            for (int i = 0; i < 50; i++)
            {
                sprintf(buffer, "D:/Codes/Boolean2016/exps/data/ref_timing/t1_%d.off", i);
                names.push_back(buffer);
            }
            test(names, expr, "D:/result10.off");
        }

        expr = "0+1";
        names.clear();
        names.push_back("D:/Codes/Boolean2016/exps/data/cmp_meshworks/buddha.off");
        names.push_back("D:/Codes/Boolean2016/exps/data/cmp_meshworks/lion.off");
        test(names, expr, "D:/result.off");

        system("pause");
        return 0;
    }

    int N = 0, nMesh;
    std::string name, output, line;
    std::stringstream buffer;

    buffer = getValidStream(configFile);
    buffer >> N;
    for (int i = 0; i < N; i++)
    {
        nMesh = 0;
        names.clear();
        buffer = getValidStream(configFile);
        buffer >> nMesh;
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

        test(names, expr, output);
    }
    configFile.close();

    system("pause");
    return 0;
}