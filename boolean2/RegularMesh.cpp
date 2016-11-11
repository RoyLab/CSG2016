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

	MemoryManager* RegularMesh::memmgr;

    RegularMesh * RegularMesh::loadFromFile(const char *fileName, uint32_t id)
    {
        OffFile file;
        readOffFile(fileName, file);

		RegularMesh* result = new RegularMesh(file, id);
        return result;
    }

    void RegularMesh::writeFile(RegularMesh & mesh, const char *fileName)
    {
        size_t nVertices = xvertices().size();
        std::vector<int> idmap(nVertices, -1);
        uint32_t vCount = 0, fCount = 0;
        std::vector<Real> vSlot;
        std::vector<int> iSlot;
        
        std::vector<MyVertex::Index> tmp, iSlotBuffer;
        cyPointT tmpVec;
        bool inverse = false;
        mesh.inverseMap.emplace_back(std::numeric_limits<uint32_t>::max(), 0);
        auto rItr = mesh.inverseMap.begin();

        for (int i = 0; i < mesh.faces().size(); i++)
        {
            if (rItr->first == i)
                inverse = true;
            else if (rItr->second == i)
            {
                ++rItr;
                if (rItr->first == i)
                    inverse = true;
                else
                    inverse = false;
            }

            IPolygon* face = mesh.faces()[i];
            if (!face->isValid()) continue;

            face->getVertices(tmp);
            iSlot.push_back(face->degree());
            for (auto idx : tmp)
            {
                uint32_t cId = idmap[idx];
                if (cId == -1)
                {
                    cId = idmap[idx] = vCount++;
                    MyVertex& v = xvertex(idx);
                    if (v.isPlaneRep())
                        tmpVec = v.ppoint().toVertexBased();
                    else
                        tmpVec = v.point();

                    XR::invCoords(mesh.m_center, mesh.m_scale, tmpVec);
                    vSlot.insert(vSlot.end(), (Real*)&tmpVec, (Real*)&tmpVec + 3);
                }
                iSlotBuffer.push_back(cId);
            }
            if (inverse)
                iSlot.insert(iSlot.end(), iSlotBuffer.rbegin(), iSlotBuffer.rend());
            else
                iSlot.insert(iSlot.end(), iSlotBuffer.begin(), iSlotBuffer.end());

            tmp.clear(); iSlotBuffer.clear();
            fCount++;
        }

        OffFile file;
        file.nEdges = 0;
        file.nFaces = fCount;
        file.nVertices = vCount;

        Real* vtmp = new Real[vSlot.size()];
        memcpy(vtmp, vSlot.data(), sizeof(Real) * vSlot.size());
        file.vertices.reset(vtmp);

        int* itmp = new int[iSlot.size()];
        memcpy(itmp, iSlot.data(), sizeof(int) * iSlot.size());
        file.indices.reset(itmp);
        
        writeOffFile(fileName, file);
    }

	RegularMesh::RegularMesh(const OffFile& file, uint32_t id)
	{
        new(this)RegularMesh();
        m_id = id;

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
			face = new Triangle(id, triId++);
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

    Oriented_side MyVertex::orientation(const XPlane & p) const
    {
        if (rep > 0)
        {
            return p.orientation(xpoint(rep - 1));
        }
        else
        {
            return p.orientation(xppoint((-rep) - 1));
        }
    }

    bool MyVertex::operator==(const XPoint & p) const
    {
        if (isPlaneRep())
            return xcppoint(id()) == p;
        else
            return p == xcpoint(id());
    }

    const XPoint & MyVertex::ppoint() const
    {
        if (!isValid() || !isPlaneRep())
            throw std::exception();

        return xppoint(id());
    }

    const cyPointT & MyVertex::point() const
    {
        if (!isValid() || isPlaneRep())
            throw std::exception();

        return xpoint(id());
    }

    MyEdge::~MyEdge()
    {
        //SAFE_DELETE(inscts);
        //SAFE_DELETE(neighbor);
    }

    void MyEdge::addAjacentFace(MyVertex::Index s, MyVertex::Index e, IPolygon * fPtr)
	{
        int ori = 1;
		if (s != ends[0]) ori = -1;

        FaceIterator itr(*this);
        for (; itr; ++itr)
        {
            if (!itr.face())
            {
                itr.faceHandle().ptr = fPtr;
                itr.faceHandle().orientation = ori;
                return;
            }
        }

        //assert(0);
        extrafhs.push_back(FH{ ori, fPtr });
	}

    int MyEdge::faceOrientation(const IPolygon * ptr) const
    {
        if (fhs[0].ptr == ptr) return fhs[0].orientation;
        if (fhs[1].ptr == ptr) return fhs[1].orientation;

        for (auto& fh : extrafhs)
            if (fh.ptr == ptr) return fh.orientation;

        return 0;
    }

    bool MyEdge::remove(const IPolygon *ptr)
    {
        FaceIterator itr(*this);
        while (itr)
        {
            if (itr.face() == ptr)
            {
                itr.faceHandle().ptr = nullptr;
                return true;
            }
            ++itr;
        }
        return false;
    }

    uint32_t MyEdge::faceCount() const
    {
        auto fItr = ConstFaceIterator(*this);
        uint32_t count = 0;
        for (; fItr; ++fItr)
        {
            if (fItr.face()->isValid())
                count++;
        }
        return count;
    }

    MyVertex & MyEdge::theOtherVertex(MyVertex::Index thiz) const
    {
        return xvertex(theOtherVId(thiz));
    }

    MyVertex::Index MyEdge::theOtherVId(MyVertex::Index thiz) const
    {
        if (ends[0] == thiz)
        {
            return ends[1];
        }
        else
        {
            assert(ends[1] == thiz);
            return ends[0];
        }
    }

    Triangle::~Triangle()
    {
        SAFE_DELETE(inscts);
    }

    cyPointT & Triangle::point(int i) const
    {
        // vIds must be non-plane rep vertex, +1, and -1, then = 0
        return xpoint(vIds[i]);
    }

    void Triangle::getVertices(std::vector<MyVertex::Index>& output) const
    {
        assert(output.empty());
        output.resize(3);
        for (int i = 0; i < 3; i++)
            output[i] = vIds[i];
    }

    bool Triangle::coherentEdge(int whichEdge) const
    {
        MyEdge& thisEdge = edge(whichEdge);
        if (thisEdge.ends[0] == vertexId((whichEdge + 1) % 3))
            return true;
        else return false;
    }

    MyVertex::Index Triangle::getTheOtherVertex(MyEdge::Index eId) const
    {
        for (int i = 0; i < 3; i++)
        {
            if (eIds[i] == eId)
                return vIds[i];
        }
        throw std::exception("cannot find edge");
    }

    uint32_t Triangle::findVertex(const XPoint & pt, PosTag tag, uint32_t *&slot)
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
            return -1;
        }
    }

    void Triangle::calcSupportingPlane()
    {
        if (!sPlane.isValid())
        {
            sPlane = XPlane(xcpoint(vIds[0]), xcpoint(vIds[1]), xcpoint(vIds[2]));
            assert(sPlane.has_on(xcpoint(vIds[0])));
            assert(sPlane.has_on(xcpoint(vIds[1])));
            assert(sPlane.has_on(xcpoint(vIds[2])));
        }
    }

    void Triangle::calcBoundingPlane()
    {
        if (!bPlanes[0].isValid())
        {
            assert(sPlane.isValid());
            const cyPointT* tmp = 
                reinterpret_cast<const cyPointT*>(sPlane.normal());
            auto normal = *tmp / tmp->Length() * 0.1;
            fp_filter_edge(reinterpret_cast<Real*>(&normal));

            auto& p = xpoint(vIds[0]);
            auto& q = xpoint(vIds[1]);
            auto& r = xpoint(vIds[2]);

            bPlanes[0].setFromPEE(q, normal, r - q);
            bPlanes[1].setFromPEE(r, normal, p - r);
            bPlanes[2].setFromPEE(p, normal, q - p);

            assert(bPlanes[0].orientation(p) == ON_POSITIVE_SIDE);
            assert(bPlanes[0].orientation(q) == ON_ORIENTED_BOUNDARY);
            assert(bPlanes[0].orientation(r) == ON_ORIENTED_BOUNDARY);
            assert(bPlanes[1].orientation(q) == ON_POSITIVE_SIDE);
            assert(bPlanes[1].orientation(p) == ON_ORIENTED_BOUNDARY);
            assert(bPlanes[1].orientation(r) == ON_ORIENTED_BOUNDARY);
            assert(bPlanes[2].orientation(r) == ON_POSITIVE_SIDE);
            assert(bPlanes[2].orientation(p) == ON_ORIENTED_BOUNDARY);
            assert(bPlanes[2].orientation(q) == ON_ORIENTED_BOUNDARY);

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

    void SubPolygon::getVertices(std::vector<MyVertex::Index>& output) const
    {
        assert(output.empty());
        output.resize(degree());
        for (uint32_t i = 0; i < degree(); i++)
            output[i] = vIds[i];
    }

    void IPolygon::getEdges(std::vector<MyEdge::Index>& output) const
    {
        assert(output.empty());
        output.resize(degree());
        for (size_t i = 0; i < degree(); i++)
            output[i] = edgeId(i);
    }

    MyEdge & IPolygon::edge(int i) const
    {
        return xedge(edgeId(i));
    }

    MyVertex& IPolygon::vertex(int i) const
    {
        return xvertex(vertexId(i));
    }

    CGALTriangle convertToCGALTriangle(const Triangle* pTri)
    {
        return CGALTriangle(convertToCGALPoint<CGALPoint>(pTri->point(0)),
            convertToCGALPoint<CGALPoint>(pTri->point(1)), convertToCGALPoint<CGALPoint>(pTri->point(2)));
    }
}

