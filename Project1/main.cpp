#include <iostream>
#include <fstream>
#include "boolean.h"

using namespace Boolean;
using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cerr << "prog [input] [log]\n";
        exit(1);
    }
    const char* input = argv[1];

    const char* log = nullptr;
    if (argc > 2)
    {
        log = argv[2];
    }

    bool ret = is_pwn(input);

    cout << "Mesh " << input << " is " << (ret ? "" : "not") << " a pwn mesh\n";

    if (log)
    {
        std::ofstream outfile(log);
        outfile << ret;
        outfile.close();
    }
    
    return 0;
}