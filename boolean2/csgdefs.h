#pragma once
#include <macros.h>
#include <linmath.h>
#include <cy/cyPoint.h>

#include <CGAL\Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Depick;

#include <CGAL\Simple_cartesian.h>
typedef CGAL::Simple_cartesian<double> Dipick;

#include <CGAL\Iso_cuboid_3.h>

#include <CGAL\Bbox_3.h>
typedef CGAL::Bbox_3 Bbox_3;

typedef uint64_t IndexPair;
typedef int8_t Indicator;

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

namespace Boolean
{
    typedef cyPoint3<double> cyPointT;

    enum PosTag
    {
        NONE = -1, INNER = 0x00,
        EDGE_0 = 0x01, EDGE_1 = 0x02, EDGE_2 = 0x04,
        VER_0 = 0x10, VER_1 = 0x20, VER_2 = 0x40,

        VER_MASK = 0xf0, EDGE_MASK = 0x0f
    };

    inline int edge_idx(PosTag tag)
    {
        switch (tag)
        {
        case EDGE_0:
            return 0;
        case EDGE_1:
            return 1;
        case EDGE_2:
            return 2;
        default:
        {
            char _errline[20];
            throw std::exception(ReportErrorString.c_str());
        }
        }
    }

    inline int vertex_idx(PosTag tag)
    {
        switch (tag)
        {
        case VER_0:
            return 0;
        case VER_1:
            return 1;
        case VER_2:
            return 2;
        default:
        {
            char _errline[20];
            throw std::exception(ReportErrorString.c_str());
        }
        }
    }

    inline PosTag edge_tag(int idx)
    {
        switch (idx)
        {
        case 0: return EDGE_0;
        case 1: return EDGE_1;
        case 2: return EDGE_2;
        default: return NONE;
        }
    }

    inline PosTag vertex_tag(int idx)
    {
        switch (idx)
        {
        case 0: return VER_0;
        case 1: return VER_1;
        case 2: return VER_2;
        default: return NONE;
        }
    }

    inline bool is_edge(PosTag tag)
    {
        return tag >= EDGE_0 && tag <= EDGE_2;
    }

    inline bool is_vertex(PosTag tag)
    {
        return tag >= VER_0 && tag <= VER_2;
    }

    enum Relation
    {
        REL_UNKNOWN = 0x0000,
        REL_INSIDE = 0x0001,
        REL_OUTSIDE = 0x0002,
        REL_SAME = 0x0004,
        REL_OPPOSITE = 0x0008,
        REL_ON_BOUNDARY = 0x000F,

        REL_NOT_AVAILABLE = -1 // (0xffffffff)
    };

    inline void inverseRelation(Relation& rel)
    {
        switch (rel)
        {
        case REL_INSIDE: rel = REL_OUTSIDE; break;
        case REL_OUTSIDE: rel = REL_INSIDE; break;
        case REL_SAME: rel = REL_OPPOSITE; break;
        case REL_OPPOSITE: rel = REL_SAME; break;
        case REL_ON_BOUNDARY: rel = REL_ON_BOUNDARY; break;
        default: ReportError("unexpected relation");
        }
    }

    inline void inverseRelation(Indicator& rel)
    {
        switch (rel)
        {
        case REL_INSIDE: rel = REL_OUTSIDE; break;
        case REL_OUTSIDE: rel = REL_INSIDE; break;
        case REL_SAME: rel = REL_OPPOSITE; break;
        case REL_OPPOSITE: rel = REL_SAME; break;
        default: ReportError("unexpected relation");
        }
    }


    enum CSGNodeOp
    {
        TYPE_UNKNOWN = 0,
        TYPE_UNION, TYPE_INTERSECT, TYPE_DIFF, TYPE_LEAF
    };

    const int MARK_BEGIN = 0xff; // 因为mark还用来在第三阶段标志有没有被访问过，所以这里让出256个数字用于这些工作
    enum Mark { SEEDED, SEEDED1, SEEDED2, VISITED, UNVISITED = MARK_BEGIN };
}
