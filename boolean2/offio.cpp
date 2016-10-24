#include "precompile.h"
#include <fstream>

#include "offio.h"

namespace XR
{
	bool OffFile::isValid() const
	{
		return nVertices > 0 && nFaces > 0 && nEdges >= 0 &&
			vertices && indices;
	}

	bool readOffFile(const char *fileName, OffFile &off)
	{
		int n = std::strlen(fileName);
		const char posfix[] = ".off";
		const int n1 = std::strlen(posfix);
		assert(n > n1);
		assert(std::strcmp(fileName + n - n1, posfix) == 0);

		std::ifstream file(fileName);
		assert(file.is_open());

		std::string header;
		file >> header;
		if (header != "OFF")
			throw std::exception("wrong format");

		file >> off.nVertices >> off.nFaces >> off.nEdges;

		const int nCoord = 3 * off.nVertices;
		Real *v = new Real[nCoord];
		for (int i = 0; i < nCoord; i++)
			file >> v[i];

		const int nFace = off.nFaces;
		const int nIdx = off.nFaces * 4;
		int *idx = new int[nIdx];
		for (int i = 0; i < nFace; i++)
		{
			file >> idx[4*i];
			assert(idx[4 * i] == 3);
			for (int j = 1; j < 4; j++)
				file >> idx[4*i+j];
		}

		assert(file.good());
		file.close();

		off.vertices.reset(v);
		off.indices.reset(idx);

		return true;
	}

}

