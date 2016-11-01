#include "precompile.h"
#include <unordered_set>
#include <unordered_map>

#include "adaptive.h"
#include "intersection.h"
#include "UndirectedGraph.hpp"
#include "RegularMesh.h"
#include "xmemory.h"

namespace Boolean
{
	namespace
	{
		enum Sign
		{
			UNKOWN = 0,
			INTERSECT_ON_LINE,
			INTERSECT_ON_POINT,
			COPLANAR,
			//LESS_THAN_INTERSECT_ON_POINT,
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
	Sign compute_intervals_isectline(Oriented_side d[3], const Triangle &pr,
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
				posA = posB = pr.boundingPlane((zeroId + 1) % 3);
				return INTERSECT_ON_POINT;
			}

			int posId = posCount[0];
			int idA = (posId + 1) % 3;
			int idB = (posId + 2) % 3;
			posA = pr.boundingPlane(idB);
			posB = pr.boundingPlane(idA);

			if (idA == zeroCount[0])
			{
				tagA = vertex_tag(idA);
				tagB = edge_tag(idA);
			}
			else
			{
				tagA = edge_tag(idB);
				tagB = vertex_tag(idB);
			}
			return INTERSECT_ON_LINE;
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

	Sign tri_tri_intersect(Triangle* t[2], TriTriInsctResult& result)
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
		Sign sign;

		t[0]->calcBoundingPlane();
		t[1]->calcBoundingPlane();

		XLine line(t[0]->supportingPlane(), t[1]->supportingPlane());

		//// 统一正方向cross(n0, n1)

		//sign = compute_intervals_isectline(da, *t[0], tagA[0], tagB[0], posA[0], posB[0]);
		// 规定cross(n1, n0)为正方向, 所以反过来传参数
		sign = compute_intervals_isectline(da, *t[0], tagB[0], tagA[0], posB[0], posA[0]);

		// 规定cross(n0, n1)为正方向
		sign = compute_intervals_isectline(db, *t[1], tagA[1], tagB[1], posA[1], posB[1]);

        line.linearOrderNoCheck(posA[1], posB[1]);

        posB[0].inverse();
        posB[1].inverse();

		assert(line.linearOrder(posA[0], posB[0]) >= 0);
		assert(line.linearOrder(posA[1], posB[1]) >= 0);

		int cmpA0B1 = line.linearOrderNoCheck(posA[0], posB[1]);
		if (cmpA0B1 < 0) return NOT_INTERSECT;
		if (cmpA0B1 == 0)
		{
			result.tagA[0] = result.tagB[0] = tagA[0];
			result.tagA[1] = result.tagB[1] = tagB[1];
			result.A = result.B = posA[0]; // or posB[1]
			return INTERSECT_ON_POINT;
		}

		int cmpA1B0 = line.linearOrderNoCheck(posA[1], posB[0]);
		if (cmpA1B0 < 0) return NOT_INTERSECT;
		if (cmpA1B0 == 0)
		{
			result.tagA[0] = result.tagB[0] = tagB[0];
			result.tagA[1] = result.tagB[1] = tagA[1];
			result.A = result.B = posA[1]; // or posB[0]
			return INTERSECT_ON_POINT;
		}

		int Acmp = line.linearOrderNoCheck(posA[0], posA[1]);
		int Bcmp = line.linearOrderNoCheck(posB[0], posB[1]);

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

        if (line.linearOrderNoCheck(result.A, result.B) == 0)
            return INTERSECT_ON_POINT;

        assert(line.linearOrderNoCheck(result.A, result.B) > 0);
		return INTERSECT_ON_LINE;
	}

	uint32_t addAndMerge(Triangle* fh[2], XPoint& pt, PosTag tag[2])
	{
		auto pMem = MemoryManager::getInstance();
		uint32_t id[2];
        uint32_t *slots[2];
		id[0] = fh[0]->findVertex(pt, tag[0], slots[0]);
		id[1] = fh[1]->findVertex(pt, tag[1], slots[1]);

        if (id[0] == id[1])
        {
            if (id[0] == INVALID_UINT32)
            {
                // add new
                uint32_t vid = pMem->insertVertex(pt);
                *(slots[0]) = *(slots[1]) = vid;
#ifdef PREP_DEBUG_INFO
                xvertex(vid).refcount += 2;
#endif
                return vid;
            }
            else return id[0];
        }
        else
        {
            uint32_t vid = std::min(id[0], id[1]);
#ifdef PREP_DEBUG_INFO
            xvertex(vid).refcount += 2;
            if (id[0] != INVALID_UINT32)
                xvertex(id[0]).refcount --;
            if (id[1] != INVALID_UINT32)
                xvertex(id[1]).refcount --;
#endif
            *(slots[0]) = *(slots[1]) = vid;
            return vid;
        }
	}

	bool insctTest(Triangle* fh0, Triangle* fh1, TriIdSet* overlaps, uint32_t meshId[2])
	{
		int id0 = fh0->id();
		int id1 = fh1->id();

		/* 统一正方向cross(n0, n1) */
		Triangle* t[2] = { fh0, fh1};
		TriTriInsctResult insctRes;
		Sign sres = tri_tri_intersect(t, insctRes);

		if (sres == NOT_INTERSECT || sres == COPLANAR)
			return false;

        fh0->addTo(intersectTriangles());
        fh1->addTo(intersectTriangles());

		uint32_t v[2];
		XPoint A(fh0->supportingPlane(), fh1->supportingPlane(), insctRes.A);
		v[0] = addAndMerge(t, A, insctRes.tagA);
		if (sres != INTERSECT_ON_POINT)
		{
			XPoint B(fh0->supportingPlane(), fh1->supportingPlane(), insctRes.B);
			v[1] = addAndMerge(t, B, insctRes.tagB);

			int eId[2] = { -1, -1 };
			for (int i = 0; i < 2; i++)
				is_same_edge(insctRes.tagA[i], insctRes.tagB[i], eId[i]);

			// seems the overlap will be removed by 
			//const MyEdge& on = t[i]->edge(eId);
			//MyEdge::ConstFaceIterator itr(on);
			//IndexPair triIdPair;
			//uint32_t triId[2] = {}
			//for (; itr; itr++)
			//{
			//		MakeIndex(triId, triIdPair);
			//	overlaps->insert(triIdPair);
			//}

			NeighborInfo ninfo;
			for (int i = 0; i < 2; i++)
			{
				int i2 = (i + 1) % 2;
				ninfo.neighborMeshId = meshId[i2];
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
					EdgePBI epbi;
                    if (sign(t[i]->supportingPlane(), t[i]->boundingPlane(eId[i]), insctRes.A) < 0)
                    {
                        epbi.pends[0] = insctRes.B.opposite();
                        epbi.pends[1] = insctRes.A.opposite();
                        epbi.ends[0] = v[1];
                        epbi.ends[1] = v[0];
                    }
                    else
                    {
                        epbi.pends[0] = insctRes.A;
                        epbi.pends[1] = insctRes.B;
                        epbi.ends[0] = v[0];
                        epbi.ends[1] = v[1];
                    }
                    epbi.neighbor.push_back(ninfo);

                    assert(sign(t[i]->supportingPlane(), t[i]->boundingPlane(eId[i]), epbi.pends[0]) > 0);
                    assert(sign(t[i]->supportingPlane(), t[i]->boundingPlane(eId[i]), epbi.pends[1]) > 0);
                    assert(XLine(t[i]->supportingPlane(), t[i]->boundingPlane(eId[i])).linearOrder(epbi.pends[0], epbi.pends[1]) > 0);

					MyEdge* edge = &t[i]->edge(eId[i]);
					if (!edge->inscts)
						edge->inscts = new InsctData<EdgePBI>;
					edge->inscts->inscts[meshId[i]].push_back(epbi);
				}
				else
				{
					FacePBI fpbi;
					fpbi.vertPlane = t[i2]->supportingPlane();
                    if (sign(t[i]->supportingPlane(), fpbi.vertPlane, insctRes.A) < 0)
                    {
                        fpbi.pends[0] = insctRes.B.opposite();
                        fpbi.pends[1] = insctRes.A.opposite();
                        fpbi.ends[0] = v[1];
                        fpbi.ends[1] = v[0];
                    }
                    else
                    {
                        fpbi.pends[0] = insctRes.A;
                        fpbi.pends[1] = insctRes.B;
                        fpbi.ends[0] = v[0];
                        fpbi.ends[1] = v[1];
                    }

					fpbi.neighbor.push_back(ninfo);

                    assert(sign(t[i]->supportingPlane(), fpbi.vertPlane, fpbi.pends[0]) > 0);
                    assert(sign(t[i]->supportingPlane(), fpbi.vertPlane, fpbi.pends[1]) > 0);
                    assert(XLine(t[i]->supportingPlane(), fpbi.vertPlane).linearOrder(fpbi.pends[0], fpbi.pends[1]) > 0);

					if (!t[i]->inscts)
						t[i]->inscts = new InsctData<FacePBI>;
					t[i]->inscts->inscts[meshId[i]].push_back(fpbi);
				}
			}

		}
		return true;
	}

	void doIntersection(std::vector<RegularMesh*>& meshes, std::vector<Octree::Node*>& intersectLeaves)
	{
		AdjacentGraph *adjGraph = nullptr;
		MeshIdTriIdMap antiOverlapMap;
		antiOverlapMap.max_load_factor(0.6f);
		adjGraph = new AdjacentGraph(meshes.size());
		auto &meshList = meshes;

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

					for (auto fh0 : *meshItr[0]->second)
					{
						for (auto fh1 : *meshItr[1]->second)
						{
							uint32_t triId[2] = { fh0->id(), fh1->id() };
							IndexPair triIdPair;
							MakeIndex(triId, triIdPair);

							int id0 = fh0->id();
							int id1 = fh1->id();

							if (antiOverlapSet->find(triIdPair) != antiOverlapSet->end())
								continue;

							antiOverlapSet->insert(triIdPair);

                            if (insctTest(fh0, fh1, antiOverlapSet, meshId))
                            {
                                assert(CGAL::do_intersect(fh0->triangle(), fh1->triangle()));
								adjGraph->setValue(meshId[0], meshId[1], true);
                            }
						}
					}
				}
			}
		}
	}
}