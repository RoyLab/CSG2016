#pragma once
#include <fstream>
#include <memory>
#include "global.h"

struct OffFile
{
	std::shared_ptr<Real> vertices;
	std::shared_ptr<size_t> indices;
	size_t nVertices, nFaces, nEdges;
};

bool readOffFile(const char*, OffFile&);
bool writeOffFile(const char*, const OffFile&);
