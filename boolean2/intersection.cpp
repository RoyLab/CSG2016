#include "precompile.h"
#include <unordered_set>
#include <unordered_map>
#include "intersection.h"
#include "Octree.h"
#include "UndirectedGraph.hpp"
#include "adaptive.h"
#include "RegularMesh.h"

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
			LESS_THAN_INTERSECT_ON_POINT,
			NOT_INTERSECT
		};
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
	Sign compute_intervals_isectline(CGAL::Oriented_side d[3], Triangle &pr, 
		PosTag& tagA, PosTag& tagB, XPlane& posA, XPlane& posB)
	{
		std::vector<int> zeroCount;
		std::vector<int> posCount;
		std::vector<int> negCount;
		for (int i = 0; i < 3; i++)
		{
			switch (d[i])
			{
			case CGAL::ON_ORIENTED_BOUNDARY:
				zeroCount.push_back(i);
				break;
			case CGAL::ON_POSITIVE_SIDE:
				posCount.push_back(i);
				break;
			case CGAL::ON_NEGATIVE_SIDE:
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

			if (d[isoVId] == CGAL::ON_NEGATIVE_SIDE)
			{
				std::swap(posA, posB);
				std::swap(tagA, tagB);
			}
			return INTERSECT_ON_LINE;
		}
		else if (zeroCount.size() == 1)
		{
			if (posCount.size() == 0)
				return INTERSECT_ON_POINT;

			int posId = posCount[0];
			int idA = (posId + 2) % 3;
			int idB = (posId + 1) % 3;
			posA = pr.boundingPlane(idA);
			posB = pr.boundingPlane(idB);

			if (idA == zeroCount[0])
			{
				tagA = edge_tag(idA);
				tagB = vertex_tag(idA);
			}
			else
			{
				tagB = edge_tag(idA);
				tagA = vertex_tag(idA);
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

			if (d[isoVId] == CGAL::ON_NEGATIVE_SIDE)
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

	void makePositive(const XPlane& p, const XPlane& q, XPlane& input)
	{
		const Real* mat[3] = { p.data(), q.data(), input.data() };
		if (GS::adaptiveDet3x3Sign(mat) < 0.0)
			input = input.opposite();
	}

	struct TriTriInsctResult
	{
		PosTag tagA[2], tagB[2];
		XPoint A, B;
	};

	Sign tri_tri_intersect(const Triangle* t[2], TriTriInsctResult& result)
	{
		typedef double FT;
		typedef cyPointT Vec3d;

		CGAL::Oriented_side db[3];
		for (size_t i = 0; i < 3; i++)
			db[i] = p[0].oriented_side(t[1].vertex(i));

		if ((db[0] == db[1]) && (db[1] == db[2]))
		{
			if (db[0] == CGAL::ON_ORIENTED_BOUNDARY) return COPLANAR;
			else return NOT_INTERSECT;
		}

		CGAL::Oriented_side da[3];
		for (size_t i = 0; i < 3; i++)
			da[i] = p[1].oriented_side(t[0].vertex(i));

		if ((da[0] == da[1]) && (da[1] == da[2]))
		{
			if (da[0] == CGAL::ON_ORIENTED_BOUNDARY) return COPLANAR;
			else return NOT_INTERSECT;
		}

		for (int i = 0; i < 2; i++)
		{
			if (!fhs[i]->data->planeRep)
				fhs[i]->data->planeRep = new PBTriangle<_R>(t[i]);
		}

		PosTag tagA[2], tagB[2];
		Plane_ext<_R> posA[2], posB[2];
		Sign sign;

		// 规定cross(n1, n0)为正方向
		sign = compute_intervals_isectline(da, t[0], *fhs[0]->data->planeRep, p[1], tagA[0], tagB[0], posA[0], posB[0]);
		if (sign == INTERSECT_ON_POINT)
			return LESS_THAN_INTERSECT_ON_POINT;

		// 规定cross(n0, n1)为正方向
		sign = compute_intervals_isectline(db, t[1], *fhs[1]->data->planeRep, p[0], tagA[1], tagB[1], posA[1], posB[1]);
		if (sign == INTERSECT_ON_POINT)
			return LESS_THAN_INTERSECT_ON_POINT;

		// 统一正方向cross(n0, n1)
		std::swap(tagA[0], tagB[0]);
		std::swap(posA[0], posB[0]);

		makePositive(p[0], p[1], posA[0]);
		makePositive(p[0], p[1], posA[1]);
		makePositive(p[0], p[1], posB[0]);
		makePositive(p[0], p[1], posB[1]);

		assert(orientation(p[0], p[1], posB[0], posA[0]) > 0);
		assert(orientation(p[0], p[1], posB[1], posA[1]) > 0);

		double cmpA0B1 = orientation(p[0], p[1], posB[1], posA[0]);
		if (cmpA0B1 < 0.) return NOT_INTERSECT;
		if (cmpA0B1 == 0.) return INTERSECT_ON_POINT;

		double cmpA1B0 = orientation(p[0], p[1], posB[0], posA[1]);
		if (cmpA1B0 < 0.) return NOT_INTERSECT;
		if (cmpA1B0 == 0.) return INTERSECT_ON_POINT;

		double Acmp = orientation(p[0], p[1], posA[0], posA[1]);
		double Bcmp = orientation(p[0], p[1], posB[0], posB[1]);

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

		// if (Acmp > 0): A0在A1右边，取A0
		result->tagA[0] = Acmp >= 0.0 ? tagA[0] : defaultTag[0];
		result->tagA[1] = Acmp <= 0.0 ? tagA[1] : defaultTag[1];
		result->A = PBPoint<_R>(p[0], p[1], (Acmp >= 0.0 ? posA[0] : posA[1]));

		// if (Bcmp > 0): B0在B1右边，取B1
		result->tagB[0] = Bcmp <= 0.0 ? tagB[0] : defaultTag[0];
		result->tagB[1] = Bcmp >= 0.0 ? tagB[1] : defaultTag[1];
		result->B = PBPoint<_R>(p[0], p[1], (Bcmp >= 0.0 ? posB[1] : posB[0]));

		assert(orientation(p[0], p[1], result->B.getPlane(2), result->A.getPlane(2)) > 0);

		return INTERSECT_ON_LINE;
	}

	bool insctTest(Triangle* fh0, Triangle* fh1, TriIdSet* overlaps, uint32_t meshId[2])
	{
		static int count = 0;
		count++;

		int id0 = fh0->id();
		int id1 = fh1->id();

		K::Triangle_3 t[2] = { fh0->data->triangle, fh1->data->triangle };
		Plane_ext<K> sp[2] = { fh0->data->sp, fh1->data->sp };
		FH fhs[2] = { fh0, fh1 };
		TriTriIsectResult<K> result;

		/* 统一正方向cross(n0, n1) */
		Sign sign = tri_tri_intersect(t, sp, &result, fhs);

		if (sign == NOT_INTERSECT || sign == INTERSECT_ON_POINT ||
			sign == COPLANAR || sign == LESS_THAN_INTERSECT_ON_POINT)
			return false;

		int addwhat[2]; // 0: nothing, 1: add vertex, 2: add line
		checkManifoldEdge(fh0, fh1, overlaps, result, addwhat, meshId);

		if (addwhat[0] == 0 && addwhat[1] == 0)
			return false;

		ItstTriangle*& it0 = fh0->data->itstTri;
		ItstTriangle*& it1 = fh1->data->itstTri;

		if (!it0) it0 = new ItstTriangle(fh0);
		if (!it1) it1 = new ItstTriangle(fh1);

		int oIdA[2] = { -1, -1 }, oIdB[2] = { -1, -1 };
		VProxyItr proxyA[2], proxyB[2];

		//std::cout << "begin compare:\n";

		getVProxy(result.A, addwhat, fhs, result.tagA, oIdA, proxyA, meshId);
		getVProxy(result.B, addwhat, fhs, result.tagB, oIdB, proxyB, meshId);

		//std::cout << "end compare:\n\n";

		for (int i = 0; i < 2; i++)
		{
			if (addwhat[i] > 0)
			{
				if (oIdA[i] == -1)
					oIdA[i] = addPoint(fhs[i], result.tagA[i], proxyA[i]);

				if (oIdB[i] == -1)
					oIdB[i] = addPoint(fhs[i], result.tagB[i], proxyB[i]);
			}
		}

		// prepare the line segments
		ItstLine line;
		if (addwhat[0] > 1)
		{
			line.pts[0].vertex = proxyA[0];
			line.pts[1].vertex = proxyB[0];
			line.plane = fhs[1];
			assert(line.check(fhs[0]->data->sp));

			fhs[0]->data->itstTri->isectLines[meshId[1]].push_back(line);
			fhs[0]->data->itstTri->meshIds.insert(meshId[1]);
		}

		if (addwhat[1] > 1)
		{
			line.pts[0].vertex = proxyB[1];
			line.pts[1].vertex = proxyA[1];
			line.plane = fhs[0];
			assert(line.check(fhs[1]->data->sp));

			fhs[1]->data->itstTri->isectLines[meshId[0]].push_back(line);
			fhs[1]->data->itstTri->meshIds.insert(meshId[0]);
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
								adjGraph->setValue(meshId[0], meshId[1], true);
						}
					}
				}
			}
		}
	}
}