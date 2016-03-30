#pragma once
#include "MyMesh.h"
#include "macroutil.h"

namespace CSG
{
    class ItstGraph
    {
        DEFINE_HANDLES;

        COMMON_PROPERTY(bool, bValid);
    public:
        ItstGraph(FH fh);
        ~ItstGraph();

    };


}