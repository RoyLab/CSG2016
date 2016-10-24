#include "precompile.h"
#include <CGAL/Point_3.h>
#include <cstring>

#include "RegularMesh.h"
#include "xmemory.h"
#include "offio.h"

namespace Boolean
{
	using namespace XR;

	namespace
	{
		template <class CGALPointT>
		CGALPointT convertToCGALPoint(const cyPointT& pt)
		{
			return CGALPointT(pt.x, pt.y, pt.z);
		}
	}
	MemoryManager* RegularMesh::memmgr;

    RegularMesh * RegularMesh::loadFromFile(const char *fileName, uint32_t id)
    {
        OffFile file;
        readOffFile(fileName, file);

		RegularMesh* result = new RegularMesh(file);
		result->m_id = id;
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
		int* indices = file.indices.get();
		Triangle* face;
		int triId = 0;
		for (int i = 0; i < file.nFaces; i++)
		{
			n = indices[ptr++];
			assert(n == 3);
			face = new Triangle(triId++);
			auto tmpIdx = indices + ptr;
			ptr += n;

			for (int j = 0; j < 3; j++)
				face->vIds[j] = tmpIdx[j]+offset;

			face->eIds[2] = memmgr->getEdgeId(tmpIdx[0], tmpIdx[1], face);
			face->eIds[0] = memmgr->getEdgeId(tmpIdx[1], tmpIdx[2], face);
			face->eIds[1] = memmgr->getEdgeId(tmpIdx[2], tmpIdx[0], face);

			m_faces[i] = face;
		}
	}

	void RegularMesh::prepareBoolean()
	{
		assert(!memmgr || memmgr == MemoryManager::getInstance());
		for (auto face: m_faces)
		{
			Triangle* pTri = reinterpret_cast<Triangle*>(face);
			pTri->cgalTri = CGALTriangle(convertToCGALPoint<CGALPoint>(memmgr->points[pTri->vIds[0]]),
				convertToCGALPoint<CGALPoint>(memmgr->points[pTri->vIds[1]]), 
				convertToCGALPoint<CGALPoint>(memmgr->points[pTri->vIds[2]]));
		}
	}

	bool MyVertex::findEdge(uint32_t other, uint32_t * result) const
	{
		for (auto& item : edges)
		{
			auto &e = xedge(item);
			if (e.ends[0] == other || e.ends[1] == other)
			{
				if (result)
					*result = item;
				return true;
			}
		}

		return false;
	}

	void MyEdge::addAjacentFace(uint32_t s, uint32_t e, IPolygon * fPtr)
	{
		FH fh;
		if (s == ends[0]) fh.orientation = 1;
		else fh.orientation = -1;

		if (!fhs[0].ptr) fhs[0] = fh;
		else if (!fhs[1].ptr) fhs[1] = fh;
		else extrafhs.push_back(fh);
	}
}

