#pragma once
#include <fstream>
#include <memory>
#include "global.h"

struct OffFile
{
	std::shared_ptr<Real> vertices;
	std::shared_ptr<size_t> indices;
	int nVertices = -1, nFaces = -1, nEdges = -1;

	bool isValid() const;
};

bool readOffFile(const char*, OffFile&);
bool writeOffFile(const char*, const OffFile&);
