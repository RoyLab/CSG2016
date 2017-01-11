#include "precompile.h"
#include <unordered_set>
#include <unordered_map>

#include <xlogger.h>
#include "adaptive.h"
#include "Octree.h"
#include "intersection.h"
#include "UndirectedGraph.hpp"
#include "RegularMesh.h"
#include "xmemory.h"

namespace Boolean
{
	namespace
	{
		enum IntersectionType
		{
			UNKOWN = 0,
			INTERSECT_ON_LINE,
			INTERSECT_ON_POINT,
			COPLANAR,
			NOT_INTERSECT
		};

		inline bool is_same_edge(PosTag tag0, PosTag tag1, int& eId)
		{
			if (tag0 == tag1)
				if (is_edge(tag0))
				{
					eId = edge_idx(tag0);
					return true;
				}

			if (is_edge(tag0) && is_vertex(tag1))
			{
				if (edge_idx(tag0) != vertex_idx(tag1))
				{
					eId = edge_idx(tag0);
					return true;
				}
				else return false;
			}
			else if (is_edge(tag1) && is_vertex(tag0))
			{
				if (edge_idx(tag1) != vertex_idx(tag0))
				{
					eId = edge_idx(tag1);
					return true;
				}
				else return false;
			}
			else if (is_vertex(tag0) && is_vertex(tag1))
			{
				int id0 = vertex_idx(tag0);
				int id1 = vertex_idx(tag1);
				eId = (id0 + 1) % 3 == id1 ? (id1 + 1) % 3 : (id0 + 1) % 3;
				return true;
			}
			else
			{
				return false;
			}
		}

	}

	class AdjacentGraph :
		public XR::UndirectedGraph<bool>
	{
	public:
		AdjacentGraph(size_t n) :XR::UndirectedGraph<bool>(n) {}
		void getIntersectPrimitives(int meshId, std::vector<int>& prims);
	};

	void AdjacentGraph::getIntersectPrimitives(int meshId, std::vector<int>& prims)
	{
		for (size_t i = 0; i < m_sz; i++)
		{
			if (getValue(meshId, i))
				prims.push_back(i);
		}
	}

	typedef std::unordered_set<IndexPair> TriIdSet;
	typedef std::unordered_map<IndexPair, TriIdSet*> MeshIdTriIdMap;

	/* A, B 的方向跟cross(ref, triangle)的方向一致 */
	IntersectionType compute_intervals_isectline(Oriented_side d[3], const Triangle &pr,
		PosTag& tagA, PosTag& tagB, XPlane& posA, XPlane& posB)
	{
		std::vector<int> zeroCount;
		std::vector<int> posCount;
		std::vector<int> negCount;
		for (int i = 0; i < 3; i++)
		{
			switch (d[i])
			{
			case ON_ORIENTED_BOUNDARY:
				zeroCount.push_back(i);
				break;
			case ON_POSITIVE_SIDE:
				posCount.push_back(i);
				break;
			case ON_NEGATIVE_SIDE:
				negCount.push_back(i);
				break;
			default:
				ReportError("");
				break;
			}
		}

		if (zeroCount.size() == 0)
		{
			int isoVId = (posCount.size() == 1 ? posCount[0] : negCount[0]);
			posA = pr.boundingPlane((isoVId + 2) % 3);
			posB = pr.boundingPlane((isoVId + 1) % 3);
			tagA = edge_tag((isoVId + 2) % 3);
			tagB = edge_tag((isoVId + 1) % 3);

			if (d[isoVId] == ON_NEGATIVE_SIDE)
			{
				std::swap(posA, posB);
				std::swap(tagA, tagB);
			}
			return INTERSECT_ON_LINE;
		}
		else if (zeroCount.size() == 1)
		{
			if (posCount.size() != 1)
			{
				int zeroId = zeroCount[0];
				tagA = tagB = vertex_tag(zeroId);
                int id1 = (zeroId + 1) % 3;
                int id2 = (zeroId + 2) % 3;

                if (posCount.size() == 2)
                {
                    posA = pr.boundingPlane(id1);
                    posB = pr.boundingPlane(id2);
                }
                else
                {
                    posA = pr.boundingPlane(id2);
                    posB = pr.boundingPlane(id1);
                }
				return INTERSECT_ON_POINT;
			}
            else
            {
                int posId = posCount[0];
                int id1 = (posId + 1) % 3;
                int id2 = (posId + 2) % 3;
                posA = pr.boundingPlane(id2);
                posB = pr.boundingPlane(id1);

                if (id1 == zeroCount[0])
                {
                    tagA = vertex_tag(id1);
                    tagB = edge_tag(id1);
                }
                else
                {
                    tagA = edge_tag(id2);
                    tagB = vertex_tag(id2);
                }
                return INTERSECT_ON_LINE;
            }
		}
		else if (zeroCount.size() == 2)
		{
			int isoVId = (posCount.size() == 1 ? posCount[0] : negCount[0]);
			posA = pr.boundingPlane((isoVId + 2) % 3);
			posB = pr.boundingPlane((isoVId + 1) % 3);

			tagA = vertex_tag((isoVId + 1) % 3);
			tagB = vertex_tag((isoVId + 2) % 3);

			if (d[isoVId] == ON_NEGATIVE_SIDE)
			{
				std::swap(posA, posB);
				std::swap(tagA, tagB);
			}
			return INTERSECT_ON_LINE;
		}
		else
		{
			ReportError("");
			return NOT_INTERSECT;
		}
	}

	struct TriTriInsctResult
	{
		PosTag tagA[2], tagB[2];
		XPlane A, B;
	};

	IntersectionType tri_tri_intersect(Triangle* t[2], TriTriInsctResult& result)
	{
		t[0]->calcSupportingPlane();
		Oriented_side db[3];
		for (size_t i = 0; i < 3; i++)
			db[i] = t[0]->supportingPlane().orientation(t[1]->point(i));

		if ((db[0] == db[1]) && (db[1] == db[2]))
		{
			if (db[0] == ON_ORIENTED_BOUNDARY) return COPLANAR;
			else return NOT_INTERSECT;
		}

		t[1]->calcSupportingPlane();
		Oriented_side da[3];
		for (size_t i = 0; i < 3; i++)
			da[i] = t[1]->supportingPlane().orientation(t[0]->point(i));

		if ((da[0] == da[1]) && (da[1] == da[2]))
		{
			if (da[0] == ON_ORIENTED_BOUNDARY) return COPLANAR;
			else return NOT_INTERSECT;
		}

		PosTag tagA[2], tagB[2];
		XPlane posA[2], posB[2];
		IntersectionType sign;

		t[0]->calcBoundingPlane();
		t[1]->calcBoundingPlane();

		PlaneLine line(t[0]->supportingPlane(), t[1]->supportingPlane());

		/// 统一正方向cross(n0, n1)
		// 规定cross(n1, n0)为正方向, 所以反过来传参数
		sign = compute_intervals_isectline(da, *t[0], tagB[0], tagA[0], posB[0], posA[0]);

		// 规定cross(n0, n1)为正方向
		sign = compute_intervals_isectline(db, *t[1], tagA[1], tagB[1], posA[1], posB[1]);

        //line.linear_order_unsafe(posA[1], posB[1]);

        posB[0].inverse();
        posB[1].inverse();

		assert(line.linear_order(posA[0], posB[0]) >= 0);
		assert(line.linear_order(posA[1], posB[1]) >= 0);

		int cmpA0B1 = line.linear_order_unsafe(posA[0], posB[1]);
		if (cmpA0B1 < 0) return NOT_INTERSECT;
		if (cmpA0B1 == 0)
		{
			result.tagA[0] = result.tagB[0] = tagA[0];
			result.tagA[1] = result.tagB[1] = tagB[1];
			result.A = result.B = posA[0]; // or posB[1]
			return INTERSECT_ON_POINT;
		}

		int cmpA1B0 = line.linear_order_unsafe(posA[1], posB[0]);
		if (cmpA1B0 < 0) return NOT_INTERSECT;
		if (cmpA1B0 == 0)
		{
			result.tagA[0] = result.tagB[0] = tagB[0];
			result.tagA[1] = result.tagB[1] = tagA[1];
			result.A = result.B = posA[1]; // or posB[0]
			return INTERSECT_ON_POINT;
		}

		int Acmp = line.linear_order_unsafe(posA[0], posA[1]);
		int Bcmp = line.linear_order_unsafe(posB[0], posB[1]);

		PosTag defaultTag[2] = { INNER, INNER };
		for (size_t i = 0; i < 2; i++)
		{
			if (is_vertex(tagA[i]) && is_vertex(tagB[i]))
			{
				int ida = vertex_idx(tagA[i]);
				int idb = vertex_idx(tagB[i]);
				int biggerId = (ida + 1) % 3 == idb ? idb : ida;
				defaultTag[i] = edge_tag((biggerId + 1) % 3);
			}
		}

		// if (Acmp < 0): A0在A1右边，取A0
		result.tagA[0] = Acmp <= 0 ? tagA[0] : defaultTag[0];
		result.tagA[1] = Acmp >= 0 ? tagA[1] : defaultTag[1];
		result.A = Acmp <= 0 ? posA[0] : posA[1];

		// if (Bcmp < 0): B0在B1右边，取B1
		result.tagB[0] = Bcmp >= 0 ? tagB[0] : defaultTag[0];
		result.tagB[1] = Bcmp <= 0 ? tagB[1] : defaultTag[1];
		result.B = Bcmp <= 0 ? posB[1] : posB[0];

        if (line.linear_order_unsafe(result.A, result.B) == 0)
            return INTERSECT_ON_POINT;

        assert(line.linear_order_unsafe(result.A, result.B) > 0);
		return INTERSECT_ON_LINE;
	}

    void add_with_hint(Triangle *fh, Triangle* fh2, const PlanePoint& pt,
        VertexIndex *&slot, PosTag tag, PosTag tag2, int i2, VertexIndex hint)
    {
        if (tag == INNER)
        {
            assert(tag2 != INNER);

            EdgeIndex eIdx = INVALID_UINT32;
            if (is_edge(tag2))
            {
                eIdx = fh2->edgeId(edge_idx(tag2));
            }
            else
            {
                int vertex_order = vertex_idx(tag2);
                if (fh2->vertex((vertex_order + 1) % 3).has_on(fh->supportingPlane()))
                {
                    eIdx = fh2->edgeId((vertex_order + 1) % 3);
                }
                else
                {
                    eIdx = fh2->edgeId((vertex_order + 2) % 3);
                }
            }

            fh->findFaceVertex(pt, eIdx,
                PlaneLine(fh2->supportingPlane(), pt.plane(2)), slot, &hint);
        }
        else
        {
            if (is_edge(tag))
            {
                if (tag2 == INNER)
                {
                    fh->findEdgeVertex(pt, tag, pt.plane(i2), slot, &hint);
                }
                else
                {
                    fh->findEdgeVertex(pt, tag, pt.plane(i2), pt.plane(2), slot, &hint);
                }
            }
            else
            {
                fh->findCornerVertex(tag, slot);
            }
        }
    }

    void add_without_hint(Triangle *fh, Triangle* fh2, const PlanePoint& pt, 
        VertexIndex *&slot, PosTag tag, PosTag tag2, int i2)
    {
        if (tag == INNER)
        {
            assert(tag2 != INNER);

            EdgeIndex eIdx = INVALID_UINT32;
            if (is_edge(tag2))
            {
                eIdx = fh2->edgeId(edge_idx(tag2));
            }
            else
            {
                int vertex_order = vertex_idx(tag2);
                if (fh2->vertex((vertex_order + 1) % 3).has_on(fh->supportingPlane()))
                {
                    eIdx = fh2->edgeId((vertex_order + 1) % 3);
                }
                else
                {
                    eIdx = fh2->edgeId((vertex_order + 2) % 3);
                }
            }

            fh->findFaceVertex(pt, eIdx,
                PlaneLine(fh2->supportingPlane(), pt.plane(2)), slot);
        }
        else
        {
            if (is_edge(tag))
            {
                if (tag2 == INNER)
                {
                    fh->findEdgeVertex(pt, tag, pt.plane(i2), slot);
                }
                else
                {
                    fh->findEdgeVertex(pt, tag, pt.plane(i2), pt.plane(2), slot);
                }
            }
            else
            {
                fh->findCornerVertex(tag, slot);
            }
        }
    }

    // the three planes of pt should be: fh's splane, fh2's splane, a boundary plane
	void addAndMerge(Triangle* fh[2], PlanePoint& pt, PosTag tag[2], VertexIndex id[2])
	{
		auto pMem = GlobalData::getObject();
        VertexIndex *slots[2] = {nullptr};

        int type[2] = { 0 };
        for (int i = 0; i < 2; ++i)
        {
            if (tag[i] == INNER)
            {
                type[i] = 3;
            }
            else if (is_edge(tag[i]))
            {
                type[i] = 2;
            }
            else
            {
                type[i] = 1;
            }
        }

        int this_id = type[0] > type[1] ? 1 : 0;
        int i2 = (this_id == 0) ? 1 : 0;

        add_without_hint(fh[this_id], fh[i2], pt, slots[this_id],
            tag[this_id], tag[i2], i2);

        id[this_id] = *slots[this_id];

        std::swap(this_id, i2);
        if (id[i2] == INVALID_UINT32)
        {
            add_without_hint(fh[this_id], fh[i2], pt, slots[this_id],
                tag[this_id], tag[i2], i2);
        }
        else
        {
            add_with_hint(fh[this_id], fh[i2], pt, slots[this_id],
                tag[this_id], tag[i2], i2, id[i2]);
        }
        id[this_id] = *slots[this_id];

        if (id[0] == id[1])
        {
            if (id[0] == INVALID_UINT32)
            {
                // add new
                VertexIndex vid = pMem->insertVertex(pt);
                *(slots[0]) = *(slots[1]) = vid;
                id[0] = id[1] = vid;
            }
            return;
        }
        else
        {
            auto minmax_pair = std::minmax(id[0], id[1]);
            if (minmax_pair.second == INVALID_UINT32)
            {
                id[0] = id[1] = minmax_pair.first;
                *slots[0] = *slots[1] = minmax_pair.first;
                return;
            }
            else
            {
                pMem->mergeVertices(id[0], id[1]);
                return;
            }
        }
	}

	bool insctTest(Triangle* fh0, Triangle* fh1, TriIdSet* overlaps, uint32_t meshId[2])
	{
		int id0 = fh0->id();
		int id1 = fh1->id();

		/* 统一正方向cross(n0, n1) */
		Triangle* t[2] = { fh0, fh1};
		TriTriInsctResult insctRes;
		IntersectionType sres = tri_tri_intersect(t, insctRes);

		if (sres == NOT_INTERSECT || sres == COPLANAR)
			return false;

		VertexIndex v[2][2]; // v[0][0], the a's first vertex index
		PlanePoint A(fh0->supportingPlane(), fh1->supportingPlane(), insctRes.A);
		addAndMerge(t, A, insctRes.tagA, v[0]);
		if (sres != INTERSECT_ON_POINT)
		{
			PlanePoint B(fh0->supportingPlane(), fh1->supportingPlane(), insctRes.B);
			addAndMerge(t, B, insctRes.tagB, v[1]);

			int eId[2] = { -1, -1 };
			for (int i = 0; i < 2; i++)
				is_same_edge(insctRes.tagA[i], insctRes.tagB[i], eId[i]);

			NeighborInfo ninfo;
			for (int i = 0; i < 2; i++)
			{
				int i2 = (i + 1) % 2;
				//ninfo.neighborMeshId = meshId[i2];
				if (eId[i2] > -1)
				{
					ninfo.type = NeighborInfo::Edge;
					ninfo.neighborEdgeId = t[i2]->edgeId(eId[i2]);
				}
				else
				{
					ninfo.type = NeighborInfo::Face;
					ninfo.pTrangle = t[i2];
				}

				if (eId[i] > -1)
				{
					EdgePbi epbi;
					MyEdge* edge = &t[i]->edge(eId[i]);

                    // let it be the same sequence as edge vertex's
                    int sequence = 1;
                    const XPlane bplane = t[i]->boundingPlane(eId[i]);
                    XPlane vertPlane = bplane;
                    if (t[i]->coherentEdge(eId[i])) vertPlane.inverse();
                    PlaneLine edgeLine(t[i]->supportingPlane(), vertPlane);
                    if (edgeLine.dot(insctRes.A) < 0)
                    {
                        //epbi.pends[0] = insctRes.B.opposite();
                        //epbi.pends[1] = insctRes.A.opposite();
                        epbi.ends[0] = v[1][i];
                        epbi.ends[1] = v[0][i];
                    }
                    else
                    {
                        //epbi.pends[0] = insctRes.A;
                        //epbi.pends[1] = insctRes.B;
                        epbi.ends[0] = v[0][i];
                        epbi.ends[1] = v[1][i];
                    }
                    //epbi.neighbor.push_back(ninfo);
                    epbi.neighbor[meshId[i2]] = ninfo;

                    //assert(edgeLine.dot(epbi.pends[0]) > 0);
                    //assert(edgeLine.dot(epbi.pends[1]) > 0);
                    //assert(edgeLine.linear_order(epbi.pends[0], epbi.pends[1]) > 0);

                    if (!edge->inscts)
                    {
                        edge->inscts = new EdgeInsctData(
                            t[i]->supportingPlane(),
                            bplane,
                            t[i], edge);
                    }

                    XR_assert(edge->inscts->checkOrientation(&epbi));
					edge->inscts->inscts[meshId[i]].push_back(epbi);
				}
				else
				{
					FacePbi fpbi;
					fpbi.vertPlane = t[i2]->supportingPlane();
                    if (orientation(t[i]->supportingPlane(), fpbi.vertPlane, insctRes.A) < 0)
                    {
                        fpbi.pends[0] = insctRes.B.opposite();
                        fpbi.pends[1] = insctRes.A.opposite();
                        fpbi.ends[0] = v[1][i];
                        fpbi.ends[1] = v[0][i];
                    }
                    else
                    {
                        fpbi.pends[0] = insctRes.A;
                        fpbi.pends[1] = insctRes.B;
                        fpbi.ends[0] = v[0][i];
                        fpbi.ends[1] = v[1][i];
                    }

                    //fpbi.neighbor.push_back(ninfo);
                    fpbi.neighbor[meshId[i2]] = ninfo;

                    //assert(orientation(t[i]->supportingPlane(), fpbi.vertPlane, fpbi.pends[0]) > 0);
                    //assert(orientation(t[i]->supportingPlane(), fpbi.vertPlane, fpbi.pends[1]) > 0);
                    //assert(PlaneLine(t[i]->supportingPlane(), fpbi.vertPlane).linear_order(fpbi.pends[0], fpbi.pends[1]) > 0);
                    //assert(linear_order(PlaneLine(t[i]->supportingPlane(), fpbi.vertPlane), fpbi.ends[0], fpbi.ends[1]) > 0);

                    if (!t[i]->inscts)
                    {
                        t[i]->inscts = new FaceInsctData(t[i]);
                    }
                    assert(t[i]->inscts->checkOrientation(&fpbi));
					t[i]->inscts->inscts[meshId[i2]].push_back(fpbi);
				}
			}
		}
		return true;
	}

    bool cgalTriTriCheck(Triangle* t0, Triangle* t1)
    {
#ifdef XR_PROFILE
        XLOG_FATAL << "assert is not turnoff!";
#endif
        CGALTriangle tr0 = convertToCGALTriangle<Depick>(t0);
        CGALTriangle tr1 = convertToCGALTriangle<Depick>(t1);

        return CGAL::do_intersect(tr0, tr1);
    }

	void doIntersection(std::vector<RegularMesh*>& meshes, std::vector<Octree::Node*>& intersectLeaves, std::vector<Triangle*>& insct_triangles)
	{
		AdjacentGraph *adjGraph = nullptr;
		MeshIdTriIdMap antiOverlapMap;
		antiOverlapMap.max_load_factor(0.6f);
		adjGraph = new AdjacentGraph(meshes.size());
		auto &meshList = meshes;

        // debug info
        int check_count = 0;
        int test_count = 0;
        int collison_count = 0;

		for (Octree::Node* leaf : intersectLeaves)
		{
			auto iEnd = leaf->triTable.cend();
			decltype(leaf->triTable.begin()) meshItr[2];
			for (meshItr[0] = leaf->triTable.begin(); meshItr[0] != iEnd; ++meshItr[0])
			{
				meshItr[1] = meshItr[0]; ++meshItr[1];
				for (; meshItr[1] != iEnd; ++meshItr[1])
				{
					uint32_t meshId[2] = { meshItr[0]->first, meshItr[1]->first };
					RegularMesh* meshes[2] = { meshList[meshId[0]], meshList[meshId[1]] };

					// 这里map的遍历是保序的，因此meshIdPair自动的分为大小
					IndexPair meshIdPair;
					MakeIndex(meshId, meshIdPair);

					TriIdSet* antiOverlapSet = nullptr;
					auto searchRes = antiOverlapMap.find(meshIdPair);
					if (searchRes == antiOverlapMap.end())
					{
						antiOverlapSet = new TriIdSet;
						antiOverlapSet->max_load_factor(0.6f);
						antiOverlapMap.emplace(meshIdPair, antiOverlapSet);
					}
					else antiOverlapSet = searchRes->second;

					for (auto fh0 : meshItr[0]->second)
					{
						for (auto fh1 : meshItr[1]->second)
						{
							uint32_t triId[2] = { fh0->id(), fh1->id() };
							IndexPair triIdPair;
							MakeIndex(triId, triIdPair);

							int id0 = fh0->id();
							int id1 = fh1->id();

                            ++check_count;
							if (antiOverlapSet->find(triIdPair) != antiOverlapSet->end())
								continue;

							antiOverlapSet->insert(triIdPair);

                            if (insctTest(fh0, fh1, antiOverlapSet, meshId))
                            {
                                ++collison_count;
                                Triangle* t[2] = {fh0, fh1};
                                for (int i = 0; i < 2; i++)
                                {
                                    if (!t[i]->add_as_insct_triangle)
                                    {
                                        t[i]->add_as_insct_triangle = true;
                                        insct_triangles.push_back(t[i]);
                                    }
                                }
                                assert(cgalTriTriCheck(fh0, fh1));
								adjGraph->setValue(meshId[0], meshId[1], true);
                            }

                            
                            if (++test_count % 500000 == 0)
                            {
                                XLOG_INFO << "Processing " << check_count << " / " 
                                    << test_count << " / " << collison_count << " triangle pairs.";
                            }
						}
					}
				}
			}
		}
	}

    VertexIndex* EdgeInsctData::find_point(const PlanePoint & p)
    {
        for (auto itr = points.begin(); itr != points.end(); itr++)
        {
            if (xvertex(itr->vertex_idx).isCoincident(p))
                return &itr->vertex_idx;
        }
        return nullptr;
    }

    VertexIndex * EdgeInsctData::find_point(const XPlane & p)
    {
        
        for (auto itr = points.begin(); itr != points.end(); itr++)
        {
            if (line.linear_coincident_no_check(p, itr->plane_rep))
            {
                return &itr->vertex_idx;
            }
        }
        return nullptr;
    }

    VertexIndex* EdgeInsctData::point(const PlanePoint & p, const XPlane& plane)
    {
        assert(line.has_on(p));
        assert(plane.has_on(p));
        assert(line.dot(plane) != 0);

        VertexIndex* result = nullptr;
        result = find_point(plane);

        if (!result)
        {
            points.push_back(Vertex{ INVALID_UINT32, plane });
            line.make_positive(points.back().plane_rep);
            result = &(points.back().vertex_idx);
        }
        return result;
    }

    auto FaceInsctData::point(const PlanePoint &p, EdgeSIndex eIdx, const PlaneLine& pline)->Vertex*
    {
        assert(pline.dot(splane_) != 0);
        for (auto itr = points.begin(); itr != points.end(); itr++)
        {
            if (orientation(splane_, pline.plane(0), itr->line.plane(0), itr->line.plane(1)) == 0)
            {
                if (orientation(splane_, pline.plane(1), itr->line.plane(0), itr->line.plane(1)) == 0)
                {
                    return &(*itr);
                }
            }
        }
        Vertex v{ INVALID_UINT32, eIdx, pline};
        assert(checkOrientation(p, eIdx));
        points.push_back(v);
        return &points.back();
    }

    bool FaceInsctData::checkOrientation(const FacePbi * pbi) const
    {
        PlaneLine line(splane_, pbi->vertPlane);
        if (linear_order(line, pbi->ends[0], pbi->ends[1]) > 0 &&
            line.dot(pbi->pends[0]) > 0 &&
            line.dot(pbi->pends[1]) > 0 &&
            line.linear_order_unsafe(pbi->pends[0], pbi->pends[1]))
        {
            return true;
        }
        else {
            return false;
        }
    }

    bool FaceInsctData::checkOrientation(const Vertex & p) const
    {
        if (!xvertex(p.vId).has_on(splane_))
        {
            return false;
        }

        if (p.eId < 0)
        {
            return true;
        }

        MyEdge& edge = xedge(p.eId);
        int a = xvertex(edge.ends[0]).has_on(splane_);
        int b = xvertex(edge.ends[1]).has_on(splane_);

        if (a + b == 2)
        {
            return false;
        }

        if (edge.inscts && !xvertex(p.vId).has_on(edge.inscts->line))
        {
            return false;
        }
        return true;
    }

    bool FaceInsctData::checkOrientation(const PlanePoint & vertex, EdgeSIndex eId) const
    {
        if (!splane_.has_on(vertex))
        {
            return false;
        }

        if (eId < 0)
        {
            return true;
        }

        MyEdge& edge = xedge(eId);
        int a = xvertex(edge.ends[0]).has_on(splane_);
        int b = xvertex(edge.ends[1]).has_on(splane_);

        if (a + b == 2)
        {
            return false;
        }

        if (edge.inscts && !edge.inscts->line.has_on(vertex))
        {
            return false;
        }
        return true;
    }

}