#include "precompile.h"
#include <CGAL/Point_3.h>
#include <cstring>

#include "RegularMesh.h"
#include "xmemory.h"
#include "offio.h"
#include "intersection.h"

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
		int n, ptr = 0;
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

			face->eIds[2] = memmgr->getEdgeId(face->vIds[0], face->vIds[1], face);
			face->eIds[0] = memmgr->getEdgeId(face->vIds[1], face->vIds[2], face);
			face->eIds[1] = memmgr->getEdgeId(face->vIds[2], face->vIds[0], face);

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

    bool MyVertex::operator==(const XPoint & p) const
    {
        if (isPlaneRep())
            return xcppoint(id()) == p;
        else
            return p == xcpoint(id());
    }

    MyEdge::~MyEdge()
    {
        SAFE_DELETE(inscts);
    }

    void MyEdge::addAjacentFace(uint32_t s, uint32_t e, IPolygon * fPtr)
	{
		FH fh( 1, fPtr );
		if (s != ends[0]) fh.orientation = -1;

		if (!fhs[0].ptr) fhs[0] = fh;
		else if (!fhs[1].ptr) fhs[1] = fh;
		else
		{
			assert(0);
			extrafhs.push_back(fh);
		}
	}

    Triangle::~Triangle()
    {
        SAFE_DELETE(inscts);
    }

    cyPointT & Triangle::point(int i) const
    {
        return xpoint(i);
    }

    MyEdge & Triangle::edge(int i) const
    {
        return xedge(eIds[i]);
    }

    int Triangle::findVertex(const XPoint & pt, PosTag tag, uint32_t *&slot)
    {
        InsctData<EdgePBI> **is;
        switch (tag)
        {
        case INNER:
            if (!inscts) inscts = new InsctData<FacePBI>;
            slot = inscts->point(pt);
            return *slot;
        case EDGE_0:
            is = &xedge(eIds[0]).inscts;
            if (!*is) *is = new InsctData<EdgePBI>;
            slot = (*is)->point(pt);
            return *slot;
        case EDGE_1:
            is = &xedge(eIds[1]).inscts;
            if (!*is) *is = new InsctData<EdgePBI>;
            slot = (*is)->point(pt);
            return *slot;
        case EDGE_2:
            is = &xedge(eIds[2]).inscts;
            if (!*is) *is = new InsctData<EdgePBI>;
            slot = (*is)->point(pt);
            return *slot;
        case VER_0:
            slot = &vIds[0];
            return *slot;
        case VER_1:
            slot = &vIds[1];
            return *slot;
        case VER_2:
            slot = &vIds[2];
            return *slot;
        default:
            assert(0);
            break;
        }
        return -1;
    }

    void Triangle::calcSupportingPlane()
    {
        if (!sPlane.isValid())
            sPlane = XPlane(xcpoint(vIds[0]), xcpoint(vIds[0]), xcpoint(vIds[0]));
    }

    void Triangle::calcBoundingPlane()
    {
        if (!bPlanes[0].isValid())
        {
            assert(sPlane.isValid());
            const cyPointT* tmp = 
                reinterpret_cast<const cyPointT*>(sPlane.normal());
            auto normal = *tmp / tmp->Length() * 0.1;
            fp_filter(reinterpret_cast<Real*>(&normal));

            auto& p = xpoint(vIds[0]);
            auto& q = xpoint(vIds[1]);
            auto& r = xpoint(vIds[2]);

            bPlanes[0] = XPlane(q, normal, r - q);
            bPlanes[1] = XPlane(r, normal, p - r);
            bPlanes[2] = XPlane(p, normal, q - p);

            assert(sPlane.has_on(p));
            assert(sPlane.has_on(q));
            assert(sPlane.has_on(r));
            assert(bPlanes[0].has_on(q));
            assert(bPlanes[0].has_on(r));
            assert(bPlanes[1].has_on(r));
            assert(bPlanes[1].has_on(p));
            assert(bPlanes[2].has_on(p));
            assert(bPlanes[2].has_on(q));
        }
    }
}

