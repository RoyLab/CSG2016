#pragma once

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
