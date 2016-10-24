#include "precompile.h"

#include <cstring>
#include "RegularMesh.h"
#include "xmemory.h"
#include "offio.h"

namespace Boolean
{
    RegularMesh * RegularMesh::loadFromFile(const char *fileName)
    {
        int n = std::strlen(fileName);
        const char posfix[] = ".off";
        const int n1 = std::strlen(posfix);
        assert(n > n1);
        assert(std::strcmp(fileName + n - n1, posfix) == 0);

        OffFile file;
        readOffFile(fileName, file);

		RegularMesh* result = new RegularMesh(file);
        return result;
    }

	RegularMesh::RegularMesh(const OffFile& file)
	{
		assert(file.isValid());
		assert(!memmgr || memmgr == MemoryManager::getInstance());
		if (!memmgr)
			memmgr = MemoryManager::getInstance();

		cyPointT* pPts = reinterpret_cast<cyPointT*>(file.vertices.get());
		int offset = memmgr->insertVertices(pPts, pPts + file.nVertices);

		m_faces.resize(file.nFaces);
		int n, ptr = 0, id;
		size_t* indices = file.indices.get();
		Triangle* face;
		for (int i = 0; i < file.nFaces; i++)
		{
			n = indices[ptr++];
			assert(n == 3);
			face = new Triangle;
			auto tmpIdx = indices + ptr;

			for (int j = 0; j < 3; j++)
				face->vIds[j] = tmpIdx[j]+offset;

			face->eIds[2] = memmgr->getEdgeId(tmpIdx[0], tmpIdx[1]);
			face->eIds[0] = memmgr->getEdgeId(tmpIdx[1], tmpIdx[2]);
			face->eIds[1] = memmgr->getEdgeId(tmpIdx[2], tmpIdx[0]);

			m_faces[i] = face;
		}
	}
}

