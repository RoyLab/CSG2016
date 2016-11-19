#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>
#include "boolean.h"

std::istringstream getValidStream(std::ifstream& file)
{
    std::string line;
    while (file)
    {
        std::getline(file, line);
        if (line.size() && line[0] != ';')
            return (std::istringstream(line));
    }
}

int main()
{
    std::string expr;
    std::vector<std::string> names;

    std::ifstream configFile(".\\mycsg.ini");
    if (!configFile.is_open())
    {
        if (false)
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
    }

    int N = 0, nMesh;
    std::string name, output, line;
    std::istringstream buffer;

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
            names.push_back(name);
        }
        buffer = getValidStream(configFile);
        buffer >> expr;
        buffer = getValidStream(configFile);
        buffer >> output;
        test(names, expr, output);
    }
    configFile.close();

    system("pause");
    return 0;
}