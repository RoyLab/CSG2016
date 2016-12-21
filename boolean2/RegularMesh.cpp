#include "precompile.h"
#include <CGAL/Point_3.h>
#include <cstring>

#include <xgeometry.h>

#include "adaptive.h"
#include "RegularMesh.h"
#include "xmemory.h"
#include "offio.h"
#include "intersection.h"

namespace Boolean
{
	using namespace XR;

	GlobalData* RegularMesh::memmgr;

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
        
        std::vector<VertexIndex> tmp, iSlotBuffer;
        cyPointT tmpVec;
        bool inverse = false;
        mesh.inverseMap.emplace_back(std::numeric_limits<uint32_t>::max(), 0);
        auto rItr = mesh.inverseMap.begin();

        for (int i = 0; i < mesh.faces().size(); i++)
        {
            if (rItr->first == i)
                inverse = true;
            else if (inverse && rItr->second == i)
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
		assert(!memmgr || memmgr == GlobalData::getObject());
		if (!memmgr)
			memmgr = GlobalData::getObject();

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
		assert(!memmgr || memmgr == GlobalData::getObject());
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

    void Triangle::getVertices(std::vector<VertexIndex>& output) const
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

    VertexIndex Triangle::getTheOtherVertex(EdgeIndex eId) const
    {
        for (int i = 0; i < 3; i++)
        {
            if (eIds[i] == eId)
                return vIds[i];
        }
        throw std::exception("cannot find edge");
    }

    VertexIndex Triangle::findVertex(const PlanePoint& pt, EdgeIndex eIdx, PosTag tag, VertexIndex*& slot)
    {
        assert(tag == INNER);
        if (!inscts) inscts = new FaceInsctData;
        slot = &inscts->point(pt, eIdx)->vId;
        return *slot;
    }


    VertexIndex Triangle::findNonFaceVertex(const PlanePoint & pt, PosTag tag, VertexIndex *&slot)
    {
        EdgeInsctData **pp_insct_data;
        switch (tag)
        {
        case EDGE_0:
            pp_insct_data = &xedge(eIds[0]).inscts;
            if (!*pp_insct_data)
            {
                *pp_insct_data = new EdgeInsctData(supportingPlane(), boundingPlane(0));
            }
            slot = (*pp_insct_data)->point(pt);
            break;
        case EDGE_1:
            pp_insct_data = &xedge(eIds[1]).inscts;
            if (!*pp_insct_data)
            {
                *pp_insct_data = new EdgeInsctData(supportingPlane(), boundingPlane(1));
            }
            slot = (*pp_insct_data)->point(pt);
            break;
        case EDGE_2:
            pp_insct_data = &xedge(eIds[2]).inscts;
            if (!*pp_insct_data) 
            {
                *pp_insct_data = new EdgeInsctData(supportingPlane(), boundingPlane(2));
            }
            slot = (*pp_insct_data)->point(pt);
            break;
        case VER_0:
            slot = &vIds[0];
            break;
        case VER_1:
            slot = &vIds[1];
            break;
        case VER_2:
            slot = &vIds[2];
            break;
        default:
            throw 1;
        }
        return *slot;

    }

    void Triangle::calcSupportingPlane()
    {
        if (!sPlane.is_valid())
        {
            sPlane = XPlane(XPlaneBase(xcpoint(vIds[0]), xcpoint(vIds[1]), xcpoint(vIds[2])));
            assert(sPlane.has_on(xcpoint(vIds[0])));
            assert(sPlane.has_on(xcpoint(vIds[1])));
            assert(sPlane.has_on(xcpoint(vIds[2])));
        }
    }

    void Triangle::calcBoundingPlane()
    {
        if (!bPlanes[0].is_valid())
        {
            assert(sPlane.is_valid());
            const cyPointT* tmp = 
                reinterpret_cast<const cyPointT*>(sPlane.normal());
            cyPointT normal = *tmp / tmp->Length() * 0.1;
            fp_filter_edge(reinterpret_cast<Real*>(&normal));

            cyPointT& p = xpoint(vIds[0]);
            cyPointT& q = xpoint(vIds[1]);
            cyPointT& r = xpoint(vIds[2]);

            bPlanes[0].setBase(XPlaneBase(q, normal, r - q, 0));
            bPlanes[1].setBase(XPlaneBase(r, normal, p - r, 0));
            bPlanes[2].setBase(XPlaneBase(p, normal, q - p, 0));

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

    void Triangle::refine()
    {
        assert(add_as_insct_triangle);

        if (inscts)
        {
            inscts->refine(this);
        }

        for (int i = 0; i < 3; i++)
        {
            if (edge(i).inscts)
            {
                edge(i).inscts->refine(this, i);
            }
        }
    }

    void SubPolygon::getVertices(std::vector<VertexIndex>& output) const
    {
        assert(output.empty());
        output.resize(degree());
        for (uint32_t i = 0; i < degree(); i++)
            output[i] = vIds[i];
    }

    void IPolygon::getEdges(std::vector<EdgeIndex>& output) const
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

