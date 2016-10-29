#include "precompile.h"
#define XRWY_EXPORTS
#include "boolean.h"
#include "preps.h"

using namespace Boolean;

extern "C"
{

    XRWY_DLL void test1()
    {
        cyPointT a(0, 0, 1), b(1, 0, 1), c(1, 1, 1);
        XPlane p(a, b, c);
        auto ptr = p.data();
        XPlane p2(a, b, c);
    }

    XRWY_DLL void test2()
    {

    }

    XRWY_DLL void test3()
    {

    }

    XRWY_DLL void test4()
    {

    }
}

