#pragma once
#include <vector>
#include "global.h"
#include "Octree.h"

namespace Boolean
{
    struct PBIRep
    {

    };

    class InsctData
    {
    public:
        typedef std::list<PBIRep> Data;

        void refine();
        void isRefined() const;
        Data& data() { return inscts; }
        const Data& data() const { return inscts; }

    protected:
        bool        bRefined = false;
        Data        inscts;
    };

    void doIntersection(std::vector<RegularMesh*>& meshes, std::vector<Octree::Node*>& intersectLeaves);
}