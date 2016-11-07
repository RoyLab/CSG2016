#pragma once
#include <cstdint>
#include "csgdefs.h"
#include "xgeometry.h"
#include "CGALext.h"

typedef double Real;

namespace Boolean
{
	typedef uint64_t IndexPair;
	typedef int8_t Indicator;

	class ICSGMesh
	{
	public:
		virtual bool& inverse() = 0;
		virtual const bool& inverse() const = 0;
		virtual uint32_t id() const = 0;
	};
    class Octree;
    class BSPTree;
	class XPlane;
	class XLine;
	class XPoint;
	class IPolygon;
	class Triangle;
    class SubPolygon;
	class RegularMesh;
	class MemoryManager;

    struct NeighborInfo;
    struct FacePBI;
    struct PBIRep;
    typedef PBIRep EdgePBI;

    template <class PBI> class InsctData;
    template <class Mesh> class CSGTree;


    const uint32_t INVALID_UINT32 = uint32_t(-1);
}
