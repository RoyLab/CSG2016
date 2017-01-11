#include "precompile.h"
#include <CGAL/Point_3.h>
#include <cstring>

#include <xgeometry.h>
#include <xlogger.h>

#include "adaptive.h"
#include "geometry.hpp"
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
        size_t nVertices = GlobalData::getObject()->vertices.size();
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

            face->get_vertices_for_dumping(tmp);
            iSlot.push_back(tmp.size());
            for (auto idx : tmp)
            {
                uint32_t cId = idmap[idx];
                if (cId == -1)
                {
                    cId = idmap[idx] = vCount++;
                    MyVertex& v = xvertex(idx);
                    if (v.isPlaneRep())
                        tmpVec = v.plane_rep().toVertexBased();
                    else
                        tmpVec = v.vertex_rep();

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

        // gen boundingbox
        bbox_ = XR::BoundingBox(pPts, pPts + file.nVertices);

        // add faces
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

			face->eIds[2] = memmgr->get_edge_id_local_scope(face->vIds[0], face->vIds[1], face);
			face->eIds[0] = memmgr->get_edge_id_local_scope(face->vIds[1], face->vIds[2], face);
			face->eIds[1] = memmgr->get_edge_id_local_scope(face->vIds[2], face->vIds[0], face);

			m_faces[i] = face;
		}
	}

    void RegularMesh::clearFaces()
    {
        for (FaceT* ptr: m_faces)
        {
            delete ptr;
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

    void Triangle::get_vertices_for_dumping(std::vector<VertexIndex>& output) const
    {
        assert(output.empty());
        output.resize(3);
        for (int i = 0; i < 3; i++)
            output[i] = vIds[i];
    }

    bool Triangle::coherentEdge(int whichEdge) const
    {
        MyEdge& thisEdge = edge(whichEdge);
        if (vertex_id_equals_simple(thisEdge.ends[0], vertexId((whichEdge + 1) % 3)))
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

    VertexIndex Triangle::get_rep_vertex(EdgeIndex id) const
    {
        MyEdge& edge = xedge(id);

        int edgeIndexInFace = -1;
        for (size_t i = 0; i < degree(); i++)
        {
            if (edgeId(i) == id)
            {
                edgeIndexInFace = i;
                break;
            }
        }
        assert(edgeIndexInFace != -1);
        return vertexId(edgeIndexInFace);
    }

    VertexIndex Triangle::findFaceVertex(const PlanePoint& pt, 
        EdgeIndex eIdx, const PlaneLine& line, VertexIndex*&slot, VertexIndex* hint)
    {
        if (!inscts)
        {
            inscts = new FaceInsctData(this);
        }
        
        bool found = false;
        if (hint)
        {
            for (auto &pts : inscts->points)
            {
                if (pts.vId == *hint)
                {
                    slot = &pts.vId;
                    found = true;
                    break;
                }
            }
        }

        if (!found)
        {
            slot = &inscts->point(pt, eIdx, line)->vId;
        }
        return *slot;
    }

    VertexIndex Triangle::findEdgeVertex(const PlanePoint & pt, PosTag tag, const XPlane& vert1,
        const XPlane& vert2, VertexIndex *&slot, VertexIndex* hint)
    {
        EdgeInsctData **pp_insct_data = nullptr;
        MyEdge* edge = nullptr;
        switch (tag)
        {
        case EDGE_0:
            edge = &xedge(eIds[0]);
            pp_insct_data = &edge->inscts;
            if (!*pp_insct_data)
            {
                *pp_insct_data = new EdgeInsctData(
                    supportingPlane(),
                    boundingPlane(0),
                    this, edge
                );
            }
            break;
        case EDGE_1:
            //pp_insct_data = &xedge(eIds[1]).inscts;
            edge = &xedge(eIds[1]);
            pp_insct_data = &edge->inscts;
            if (!*pp_insct_data)
            {
                *pp_insct_data = new EdgeInsctData(
                    supportingPlane(),
                    boundingPlane(1),
                    this, edge
                );
            }
            break;
        case EDGE_2:
            //pp_insct_data = &xedge(eIds[2]).inscts;
            edge = &xedge(eIds[2]);
            pp_insct_data = &edge->inscts;
            if (!*pp_insct_data)
            {
                *pp_insct_data = new EdgeInsctData(
                    supportingPlane(),
                    boundingPlane(2),
                    this, edge
                );
            }
            break;
        default:
            throw 1;
        }

        bool found = false;
        if (hint)
        {
            for (auto &pts : edge->inscts->points)
            {
                if (pts.vertex_idx == *hint)
                {
                    slot = &pts.vertex_idx;
                    found = true;
                    break;
                }
            }
        }

        if (!found)
        {
            if ((*pp_insct_data)->line.dot(vert1) == 0)
            {
                slot = (*pp_insct_data)->point(pt, vert2);
            }
            else
            {
                slot = (*pp_insct_data)->point(pt, vert1);
            }
        }
        return *slot;
    }

    VertexIndex Triangle::findEdgeVertex(const PlanePoint & pt, PosTag tag, 
        const XPlane& vert, VertexIndex *&slot, VertexIndex* hint)
    {
        EdgeInsctData **pp_insct_data = nullptr;
        MyEdge* edge = nullptr;
        switch (tag)
        {
        case EDGE_0:
            edge = &xedge(eIds[0]);
            pp_insct_data = &edge->inscts;
            if (!*pp_insct_data)
            {
                *pp_insct_data = new EdgeInsctData(
                    supportingPlane(),
                    boundingPlane(0),
                    this, edge
                );
            }
            break;
        case EDGE_1:
            //pp_insct_data = &xedge(eIds[1]).inscts;
            edge = &xedge(eIds[1]);
            pp_insct_data = &edge->inscts;
            if (!*pp_insct_data)
            {
                *pp_insct_data = new EdgeInsctData(
                    supportingPlane(),
                    boundingPlane(1),
                    this, edge
                );
            }
            break;
        case EDGE_2:
            //pp_insct_data = &xedge(eIds[2]).inscts;
            edge = &xedge(eIds[2]);
            pp_insct_data = &edge->inscts;
            if (!*pp_insct_data)
            {
                *pp_insct_data = new EdgeInsctData(
                    supportingPlane(),
                    boundingPlane(2),
                    this, edge
                );
            }
            break;
        default:
            throw 1;
        }

        bool found = false;
        if (hint)
        {
            for (auto &pts : edge->inscts->points)
            {
                if (pts.vertex_idx == *hint)
                {
                    slot = &pts.vertex_idx;
                    found = true;
                    break;
                }
            }
        }

        if (!found)
        {
            slot = (*pp_insct_data)->point(pt, vert);
        }
        return *slot;
    }

    VertexIndex Triangle::findCornerVertex(PosTag tag, VertexIndex *&slot)
    {
        switch(tag)
        {
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
            sPlane.construct_from_three_vertices(xcpoint(vIds[0]), xcpoint(vIds[1]), xcpoint(vIds[2]));
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
            //const cyPointT* tmp = 
            //    reinterpret_cast<const cyPointT*>(sPlane.normal());
            //cyPointT normal = *tmp / tmp->Length() * 0.1;
            //fp_filter_edge(reinterpret_cast<Real*>(&normal));
            cyPointT normal = get_projection_vector_unit(sPlane.normal()) / std::pow(2, 10);

            assert(fp_filter_check(
                reinterpret_cast<Real*>(&normal), 
                FP_EDGE_CHECK
            ));

            cyPointT& p = xpoint(vIds[0]);
            cyPointT& q = xpoint(vIds[1]);
            cyPointT& r = xpoint(vIds[2]);

            bPlanes[0].construct_from_one_vertex_two_edges(q, normal, r - q);
            bPlanes[1].construct_from_one_vertex_two_edges(r, normal, p - r);
            bPlanes[2].construct_from_one_vertex_two_edges(p, normal, q - p);

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

    void Triangle::getAllEdges(std::vector<EdgeIndex>& output) const
    {
        assert(output.empty());
        output.resize(degree());
        for (size_t i = 0; i < degree(); i++)
            output[i] = edgeId(i);
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

    void Triangle::load_coords(std::deque<Real[9]>& container)
    {
        tmp_id = container.size();
        container.emplace_back();
        cyPointT pts=point(0);
        memcpy(container.back(), &pts, sizeof(Real) * 3);
        pts = point(1);
        memcpy(container.back()+3, &pts, sizeof(Real) * 3);
        pts = point(2);
        memcpy(container.back()+6, &pts, sizeof(Real) * 3);
    }

    void SubPolygon::get_vertices_for_dumping(std::vector<VertexIndex>& output) const
    {
        assert(output.empty());
        output.resize(degree());
        for (uint32_t i = 0; i < degree(); i++)
            output[i] = vIds[i];
    }

    VertexIndex SubPolygon::get_rep_vertex(EdgeIndex id) const
    {
        assert(orientation(supportingPlane(), (vertex(0))) == ON_ORIENTED_BOUNDARY);
        assert(orientation(supportingPlane(), (vertex(1))) == ON_ORIENTED_BOUNDARY);
        assert(orientation(supportingPlane(), (vertex(2))) == ON_ORIENTED_BOUNDARY);

        int edgeIndexInFace = -1;
        for (int i = 0; i < degree(); i++)
        {
            if (edgeId(i) == id)
            {
                edgeIndexInFace = i;
                break;
            }
        }
        assert(edgeIndexInFace != -1);

        VertexIndex vIdInPlane;
        XPlane boundPlane, tmpPlane;
        MyEdge& edge = xedge(id);

        // find init point
        for (int i = 2; i < degree(); i++)
        {
            vIdInPlane = vertexId((edgeIndexInFace + i) % degree());

            // find a bounding plane
            bool flag = false;
            for (auto &neigh : *edge.neighbor)
            {
                if (neigh.second.type == NeighborInfo::Edge)
                {
                    for (auto fItr = MyEdge::FaceIterator(xedge(neigh.second.neighborEdgeId), true);
                        fItr; fItr.incrementToTriangle())
                    {
                        assert(fItr.face()->getType() == IPolygon::TRIANGLE);
                        tmpPlane = ((Triangle*)fItr.face())->supportingPlane();
                        if (orientation(tmpPlane, xvertex(vIdInPlane)) != ON_ORIENTED_BOUNDARY)
                        {
                            boundPlane = tmpPlane;
                            flag = true;
                            break;
                        }
                    }
                }
                else
                {
                    assert(neigh.second.type == NeighborInfo::Face);
                    XPlane tmpPlane = neigh.second.pTrangle->supportingPlane();
                    if (orientation(tmpPlane, xvertex(vIdInPlane)) != ON_ORIENTED_BOUNDARY)
                    {
                        boundPlane = tmpPlane;
                        break;
                    }
                }
                if (flag) break;
            }
            if (boundPlane.is_valid()) break;
        }
        assert(boundPlane.is_valid());

        // correct the direction of bounding plane
        PlaneLine edgeLine(supportingPlane(), boundPlane);
        assert(!supportingPlane().id_equals(boundPlane));
        int tmpSide = linear_order(edgeLine, xvertex(vertexId((edgeIndexInFace + 1) % degree())),
            xvertex(vertexId(edgeIndexInFace)));

        assert(tmpSide != 0);
        if (tmpSide < 0)
            boundPlane.inverse();

        // pick a correct rep vertex
        VertexIndex repVertexId = INVALID_UINT32;
        for (size_t i = 0; i < degree(); i++)
        {
            if (orientation(boundPlane, xvertex(vertexId(i))) == ON_POSITIVE_SIDE)
            {
                repVertexId = vertexId(i);
                break;
            }
        }
        if (repVertexId == INVALID_UINT32)
        {
            XLOG_ERROR << "may be there is self-intersection.";
            throw 1;
        }
        return repVertexId;
    }

    void SubPolygon::getAllEdges(std::vector<EdgeIndex>& output) const
    {
        assert(output.empty());
        output.resize(degree());
        for (size_t i = 0; i < degree(); i++)
            output[i] = edgeId(i);
    }

    MyEdge & SubPolygon::edge(int i) const
    {
        return xedge(edgeId(i));
    }

    MyVertex& SubPolygon::vertex(int i) const
    {
        return xvertex(vertexId(i));
    }

    MyEdge & Triangle::edge(int i) const
    {
        return xedge(edgeId(i));
    }

    MyVertex& Triangle::vertex(int i) const
    {
        return xvertex(vertexId(i));
    }

    SubPolygonWithHoles::SubPolygonWithHoles(const Triangle* tri, std::vector<std::vector<VertexIndex>>& loops, uint32_t i):
        IPolygon(i, tri->meshId()),  father_(tri)
    {
        loops_.resize(loops.size());
        int count = 0;
        for (auto& loop : loops)
        {
            auto itr = loop.begin();
            //VertexIndex v0, v1;
            auto pMem = GlobalData::getObject();

            loops_[count].vIds.resize(loop.size());
            loops_[count].eIds.resize(loop.size());
            for (int i = 0; i < loop.size(); i++)
            {
                //v0 = *itr; ++itr;
                //if (itr != loop.end())
                //{
                //    v1 = *itr;
                //}
                //else
                //{
                //    v1 = loop.front();
                //}

                loops_[count].vIds[i] = pMem->get_main_vertexId(*itr);
                ++itr;
            }

            for (int i = 0; i < loop.size(); ++i)
            {
                loops_[count].eIds[i] = pMem->get_edge_id_local_scope(
                    loops_[count].vIds[i], 
                    loops_[count].vIds[(i+1)% loop.size()],
                    this
                );
            }

            ++count;
        }
    }

    void SubPolygonWithHoles::get_vertices_for_dumping(std::vector<VertexIndex>& output) const
    {
        assert(output.empty());
        output = loops_[0].vIds;
    }

    void SubPolygonWithHoles::getAllEdges(std::vector<EdgeIndex>& output) const
    {
        assert(output.empty());
        for (const Loop& loop: loops_)
        { 
            output.insert(output.end(), loop.eIds.begin(), loop.eIds.end());
        }
    }

    VertexIndex SubPolygonWithHoles::get_rep_vertex(EdgeIndex id) const
    {
        assert(orientation(supportingPlane(), (vertex(0, 0))) == ON_ORIENTED_BOUNDARY);
        assert(orientation(supportingPlane(), (vertex(0, 1))) == ON_ORIENTED_BOUNDARY);
        assert(orientation(supportingPlane(), (vertex(0, 2))) == ON_ORIENTED_BOUNDARY);

        int edgeIndexInFace = -1, edgeIndexInFace2 = -1;
        for (int i = 0; i < loops_.size(); ++i)
        {
            for (int j = 0; j < loops_[i].vIds.size(); ++j)
            {
                if (edgeId(i, j) == id)
                {
                    edgeIndexInFace = i;
                    edgeIndexInFace2 = j;
                    break;
                }
            }
        }
        assert(edgeIndexInFace != -1);

        VertexIndex vIdInPlane;
        XPlane boundPlane, tmpPlane;
        MyEdge& edge = xedge(id);

        // find init point
        for (int i = 0; i < loops_.size(); ++i)
        {
            for (int j = 0; j < loops_[i].vIds.size(); ++j)
            {
                vIdInPlane = vertexId(i, j);

                // find a bounding plane
                bool flag = false;
                for (auto &neigh : *edge.neighbor)
                {
                    if (neigh.second.type == NeighborInfo::Edge)
                    {
                        for (auto fItr = MyEdge::FaceIterator(xedge(neigh.second.neighborEdgeId), true);
                            fItr; fItr.incrementToTriangle())
                        {
                            assert(fItr.face()->getType() == IPolygon::TRIANGLE);
                            tmpPlane = ((Triangle*)fItr.face())->supportingPlane();
                            if (orientation(tmpPlane, xvertex(vIdInPlane)) != ON_ORIENTED_BOUNDARY)
                            {
                                boundPlane = tmpPlane;
                                flag = true;
                                break;
                            }
                        }
                    }
                    else
                    {
                        assert(neigh.second.type == NeighborInfo::Face);
                        XPlane tmpPlane = neigh.second.pTrangle->supportingPlane();
                        if (orientation(tmpPlane, xvertex(vIdInPlane)) != ON_ORIENTED_BOUNDARY)
                        {
                            boundPlane = tmpPlane;
                            break;
                        }
                    }
                    if (flag) break;
                }
                if (boundPlane.is_valid()) break;
            }
        }
        assert(boundPlane.is_valid());

        // correct the direction of bounding plane
        PlaneLine edgeLine(supportingPlane(), boundPlane);
        assert(!supportingPlane().id_equals(boundPlane));
        int tmpSide = linear_order(edgeLine, 
            xvertex(vertexId(edgeIndexInFace, (edgeIndexInFace2 + 1) % loops_[edgeIndexInFace].vIds.size())),
            xvertex(vertexId(edgeIndexInFace, edgeIndexInFace2))
        );

        assert(tmpSide != 0);
        if (tmpSide < 0)
            boundPlane.inverse();

        // pick a correct rep vertex
        VertexIndex repVertexId = INVALID_UINT32;
        for (int i = 0; i < loops_.size(); ++i)
        {
            for (int j = 0; j < loops_[i].vIds.size(); ++j)
            {
                if (orientation(boundPlane, xvertex(vertexId(i, j))) == ON_POSITIVE_SIDE)
                {
                    repVertexId = vertexId(i, j);
                    break;
                }
            }
        }
        assert(repVertexId != INVALID_UINT32);
        return repVertexId;
    }

    MyEdge & SubPolygonWithHoles::edge(int i, int j) const
    {
        return xedge(edgeId(i, j));
    }

    MyVertex & SubPolygonWithHoles::vertex(int i, int j) const
    {
        return xvertex(vertexId(i, j));
    }
}

