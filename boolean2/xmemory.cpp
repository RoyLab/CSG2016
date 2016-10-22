#include "xmemory.h"

int Boolean::MemoryManager::insertVertices(Point * begin, Point * end)
{
	int offset = points.size();
	int n = end - begin;
	points.insert(points.end(), begin, end);

	std::vector<MyVertex> tmp(n);
	vertices.insert(vertices.end(), tmp.begin(), tmp.end());

	return offset;
}
