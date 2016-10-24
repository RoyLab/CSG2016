#include "precompile.h"
#include "xmemory.h"
#include "RegularMesh.h"


namespace Boolean
{
	int MemoryManager::insertVertices(Point * begin, Point * end)
	{
		int offset = points.size();
		int n = end - begin;

		points.insert(points.end(), begin, end);
		vertices.resize(vertices.size() + n);

		return offset;
	}

	int MemoryManager::getEdgeId(int a, int b, IPolygon * facePtr)
	{
		MyEdge edge;
		edge.addAjacentFace(a, b);
		assert(0);

		return 0;
	}
}
