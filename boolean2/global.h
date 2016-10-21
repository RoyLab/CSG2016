#pragma once

typedef double Real;

namespace Boolean
{
	enum Relation
	{
		REL_UNKNOWN = 0x0000,
		REL_INSIDE = 0x0001,
		REL_OUTSIDE = 0x0002,
		REL_SAME = 0x0004,
		REL_OPPOSITE = 0x0008,
		REL_ON_BOUNDARY = 0x000F,

		REL_NOT_AVAILABLE = -1
	};

	enum CSGNodeOp
	{
		TYPE_UNKNOWN = 0,
		TYPE_UNION, TYPE_INTERSECT, TYPE_DIFF, TYPE_LEAF
	};

	typedef uint64_t IndexPair;
	typedef int8_t Indicator;

	class XPlane;
	class XPoint;
	class RegularMesh;

	class MemoryManager;
}
