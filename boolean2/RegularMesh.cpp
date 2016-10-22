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

		MemoryManager* pMemmgr = MemoryManager::getInstance();
		cyPointT* pPts = reinterpret_cast<cyPointT*>(file.vertices.get());
		int offset = pMemmgr->insertVertices(pPts, pPts + file.nVertices);

		RegularMesh* result = new RegularMesh;

        return result;
    }
}

