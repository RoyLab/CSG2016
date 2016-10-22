#pragma once
#include <cstdint>
#include "csgdefs.h"

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
		virtual size_t id() const = 0;
	};
    class Octree;
	class XPlane;
	class XPoint;
	class RegularMesh;
	class MemoryManager;
}
