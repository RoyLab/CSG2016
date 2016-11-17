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
        std::vector<int> idx, tmpIds;
        int iPtr = 0;
        int nFacePlus = 0;
		for (int i = 0; i < nFace; i++)
		{
            int nD, tmpIdx;
            file >> nD;
            idx.push_back(3);
            tmpIds.resize(nD);

            for (int j = 0; j < nD; j++)
            {
                file >> tmpIdx;
                tmpIds[j] = tmpIdx;
            }

            switch (nD)
            {
            case 3:
                idx.insert(idx.end(), tmpIds.begin(), tmpIds.end());
                break;
            case 4:
                idx.insert(idx.end(), tmpIds.begin(), tmpIds.begin() + 3);
                idx.push_back(3);
                idx.push_back(tmpIds[0]);
                idx.push_back(tmpIds[2]);
                idx.push_back(tmpIds[3]);
                ++nFacePlus;
                break;
            default:
                throw std::exception();
            }
            tmpIds.clear();
		}

		assert(file.good());
		file.close();

        off.nFaces += nFacePlus;
		off.vertices.reset(v);

        int *idxPtr = new int[idx.size()];
        std::copy(idx.begin(), idx.end(), idxPtr);
		off.indices.reset(idxPtr);

		return true;
	}

    bool writeOffFile(const char * fileName, const OffFile &mesh)
    {
        std::ofstream f(fileName);
        if (!f.is_open()) throw std::exception("no file");

        f << "OFF\n";
        f << mesh.nVertices << ' ' << mesh.nFaces << ' ' << mesh.nEdges << std::endl;

        for (uint32_t i = 0; i < mesh.nVertices; i++)
            f << mesh.vertices.get()[3*i] << ' ' << mesh.vertices.get()[3*i+1]<<
            ' ' << mesh.vertices.get()[3*i+2] << std::endl;

        int ptr = 0;
        for (uint32_t i = 0; i < mesh.nFaces; i++)
        {
            int n = mesh.indices.get()[ptr++];
            f << n;
            for (int j = 0; j < n; j++)
                f << ' ' << mesh.indices.get()[ptr++];
            f << std::endl;
        }

        f.close();
        return false;
    }

}

