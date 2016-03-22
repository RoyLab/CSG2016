#pragma once
#include "linmath.h"

#include <CGAL\Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

#include <CGAL\Iso_cuboid_3.h>
typedef CGAL::Iso_cuboid_3<K> Cube_3;

#include <CGAL\Bbox_3.h>
typedef CGAL::Bbox_3 Bbox_3;

typedef uint64_t IndexPair;

static inline void MakeIndex(const uint32_t id[], IndexPair& indexPair)
{
    indexPair = id[1];
    indexPair = indexPair << 32;
    indexPair |= id[0];
}

static inline void GetIDFromIndex(uint32_t ID[], const IndexPair& indexPair)
{
    ID[0] = indexPair & 0xffffffff;
    ID[1] = indexPair >> 32;
}

typedef double double2[2];
typedef double double3[3];
typedef double double4[4];

typedef double2 double2x2[2];
typedef double3 double3x3[3];
typedef double4 double4x4[4];

namespace CSG
{
    enum PosTag
    {
        NONE = -1, INNER = 0x00,
        EDGE_0 = 0x01, EDGE_1 = 0x02, EDGE_2 = 0x04,
        VER_0 = 0x10, VER_1 = 0x20, VER_2 = 0x40,

        VER_MASK = 0xf0, EDGE_MASK = 0x0f
    };

    enum Relation
    {
        REL_UNKNOWN = 0x0000,
        REL_INSIDE = 0x0001,
        REL_OUTSIDE = 0x0002,
        REL_SAME = 0x0004,
        REL_OPPOSITE = 0x0008,

        REL_NOT_AVAILABLE = -1 // (0xffffffff)
    };

    inline void inverseRelation(Relation& rel)
    {
        switch (rel)
        {
        case CSG::REL_INSIDE: rel = REL_OUTSIDE;
        case CSG::REL_OUTSIDE: rel = REL_INSIDE;
        case CSG::REL_SAME: rel = REL_OPPOSITE;
        case CSG::REL_OPPOSITE: rel = REL_SAME;
        }
    }


    enum CSGNodeOp
    {
        TYPE_UNKNOWN = 0,
        TYPE_UNION, TYPE_INTERSECT, TYPE_DIFF, TYPE_LEAF
    };

}
