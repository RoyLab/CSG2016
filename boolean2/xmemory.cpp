#include "xmemory.h"

int Boolean::MemoryManager::insertVertices(Point * begin, Point * end)
{
	int offset = points.size();
	int n = end - begin;
	
	points.insert(points.end(), begin, end);
	vertices.resize(vertices.size() + n);

	return offset;
}
