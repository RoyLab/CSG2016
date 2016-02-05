#include "precompile.h"
#include "BSP2D.h"
#include "Bool.h"
#include "BaseMesh.h"
#include <queue>

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

#include "MPMesh.h"
#include <ctime>
#include "configure.h"
#include "COctree.h"
#include "IsectTriangle.h"
#include "isect.h"
#include <Fade_2D.h>

namespace CSG
{

extern HANDLE _output;
extern clock_t t0;


void DebugInfo(char* str, clock_t& t0);
void StdOutput(const char* str);


class FeitoISectZone
{
public:
    struct ISVertexInfo;
	typedef std::list<ISVertexInfo>::iterator	ISVertexItr;

	struct ISVertexInfo
	{
        Vec3d coord;
        MPMesh2::VertexHandle vhandle;
		int Id_2d;
    };

    struct SegData
    {
        MPMesh2::FaceHandle thatFace;
        ISVertexItr points[2];

		Line2D lineCoef;
    };

    struct VertexData
    {
        Vec3d coord;
        MPMesh::VertexHandle used;
        VertexData(){used.reset();}
        VertexData(const Vec3d& v):coord(v){used.reset();}
    };

    FeitoISectZone(MPMesh2* mesh, MPMesh2* mesh2, MPMesh2::FaceHandle fh):
    face(fh), thisMesh(mesh), thatMesh(mesh2), xi(-1), yi(-1)
    {
		// 三个角点在vertices的最后三个
		auto fvItr = thisMesh->fv_begin(face);
		ISVertexInfo info;
		info.vhandle = *fvItr;
		info.coord = thisMesh->point(info.vhandle);
		vertices.push_back(info);
		corner[0] = vertices.end();
		corner[0] --;

		fvItr++;
		info.vhandle = *fvItr;
		info.coord = thisMesh->point(info.vhandle);
		vertices.push_back(info);
		corner[1] = vertices.end();
		corner[1] --;

		fvItr++;
		info.vhandle = *fvItr;
		info.coord = thisMesh->point(info.vhandle);
		vertices.push_back(info);
		corner[2] = vertices.end();
		corner[2] --;

		zone_record.push_back(this);
    }

    ~FeitoISectZone(){}

    void InsertSegment(ISVertexItr& v0, ISVertexItr& v1, FeitoISectZone* tri2)
	{
		segs.emplace_back();
        auto &seg = segs.back();
        seg.points[0] = v0;
        seg.points[1] = v1;
		seg.thatFace = tri2->face;
	}

	std::list<MPMesh2::FaceHandle> coplanarTris;
	std::list<SegData> segs;
    std::list<ISVertexInfo> vertices;
    MPMesh2::FaceHandle face;

    ISVertexItr corner[3];
    MPMesh2 *thisMesh, *thatMesh;

	int xi, yi;

private:
    int Id;

public:
    static void clearList(){verticesList.clear();}

    static void clear()
    {
        for (auto itr:zone_record) SAFE_RELEASE(itr);
        zone_record.clear();
        clearList();
    }

    static VertexData& point(unsigned id){return verticesList[id];}

    static unsigned add_point(Vec3d& vec)
    {
        verticesList.emplace_back(vec);
        return verticesList.size()-1;
        return 0;
    }


    static std::vector<FeitoISectZone*> zone_record;

private:
    static std::vector<VertexData> verticesList;

};

std::vector<FeitoISectZone::VertexData> FeitoISectZone::verticesList;
std::vector<FeitoISectZone*> FeitoISectZone::zone_record;

inline void GetCorners(MPMesh2* pMesh, MPMesh2::FaceHandle fhandle, Vec3d*& v0, Vec3d *&v1, Vec3d *&v2)
{
	auto fvItr = pMesh->fv_begin(fhandle);
	v0 = &(pMesh->verticesList[fvItr->idx()]);fvItr++;
	v1 = &(pMesh->verticesList[fvItr->idx()]);fvItr++;
	v2 = &(pMesh->verticesList[fvItr->idx()]);
}

void GetLeafNodes(OctreeNode*, std::list<OctreeNode*>&, int);

inline bool CompareVertex(FeitoISectZone::ISVertexItr& ref, const Vec3d& vec)
{
	return IsEqual(ref->coord, vec);
}

bool IsVertexExisted(FeitoISectZone* tri, Vec3d& vec, FeitoISectZone::ISVertexItr& ref)
{
	const auto end = tri->corner[0];
	for (auto itr = tri->vertices.begin(); itr != end; itr++)
	{
		if (CompareVertex(itr, vec))
		{
			ref = itr;
			return true;
		}
	}
	return false;
}

// 明确一点：ref传进来之后将取代同位置的点
void InsertPoint(FeitoISectZone* tri, VertexPos pos, FeitoISectZone::ISVertexInfo& ref)
{
    Vec3d &vec = ref.coord;
    FeitoISectZone::ISVertexItr coincident;
	if (pos == INNER)
	{
		if (!IsVertexExisted(tri, vec, coincident))
			tri->vertices.push_front(ref);
		else
            coincident->vhandle = ref.vhandle;
	}
	else if (pos <= EDGE_2)
	{
		auto mesh = tri->thisMesh;
		auto ffItr = mesh->ff_begin(tri->face);
		switch (pos)
		{
		case CSG::EDGE_0:
			ffItr ++;
			ffItr ++;
			break;
		case CSG::EDGE_1:
			break;
		case CSG::EDGE_2:
			ffItr ++;
			break;
		default:
			assert(0);
		}

		if (ffItr->is_valid())
		{
			FeitoISectZone*& other = mesh->property(mesh->SurfacePropHandle, *ffItr);
			if (!other) other = new FeitoISectZone(mesh, tri->thatMesh, *ffItr);
			InsertPoint(other, INNER, ref);
		}
		InsertPoint(tri, INNER, ref);
	}
	else
	{
		// get the iterator of vertex #WR#
        assert(0);
		    //switch (pos)
		    //{
		    //case CSG::VER_0:
			   // if (tri->corner[0]->pos == ref->pos) break;
			   // tri->corner[0]->reset();
			   // tri->corner[0]->next = ref;
			   // output = tri->corner[0];
			   // break;
		    //case CSG::VER_1:
			   // if (tri->corner[1]->pos == ref->pos) break;
			   // tri->corner[1]->reset();
			   // tri->corner[1]->next = ref;
			   // output = tri->corner[1];
			   // break;
		    //case CSG::VER_2:
			   // if (tri->corner[2]->pos == ref->pos) break;
			   // tri->corner[2]->reset();
			   // tri->corner[2]->next = ref;
			   // output = tri->corner[2];
			   // break;
		    //}
	}
}

FeitoISectZone::ISVertexItr InsertPoint(FeitoISectZone* tri, VertexPos pos, Vec3d& vec)
{
	FeitoISectZone::ISVertexItr output;
    FeitoISectZone::ISVertexInfo info;
    info.coord = vec;
	if (pos == INNER)
	{
		if (!IsVertexExisted(tri, vec, output))
		{
            auto vhandle = tri->thisMesh->add_vertex(vec);
            info.vhandle = vhandle;
			tri->vertices.push_front(info);
            return tri->vertices.begin();
		}
        else return output;
	}
	else if (pos <= EDGE_2)
	{
		auto mesh = tri->thisMesh;
		auto ffItr = mesh->ff_begin(tri->face);
		switch (pos)
		{
		case CSG::EDGE_0:
			ffItr ++;
			ffItr ++;
			break;
		case CSG::EDGE_1:
			break;
		case CSG::EDGE_2:
			ffItr ++;
			break;
		default: 
			assert(0);
		}

		output = InsertPoint(tri, INNER, vec);
		if (ffItr->is_valid())
		{
			FeitoISectZone*& other = mesh->property(mesh->SurfacePropHandle, *ffItr);
#ifdef _DEBUG_
			auto fvItr = mesh->fv_begin(tri->face);
			if (pos == EDGE_1) fvItr++;
			if (pos == EDGE_2) {fvItr++; fvItr++;}
			
			auto f0 = *fvItr;
			fvItr = mesh->fv_begin(*ffItr);
			assert(*fvItr != f0); fvItr++;
			assert(*fvItr != f0); fvItr++;
			assert(*fvItr != f0);
#endif
			if (!other) other = new FeitoISectZone(mesh, tri->thatMesh, *ffItr);
			InsertPoint(other, INNER, *output);
		}
        return output;
	}
	else
	{
		// get the iterator of vertex #WR#
		switch (pos)
		{
		case CSG::VER_0:
			output = tri->corner[0];
			break;
		case CSG::VER_1:
			output = tri->corner[1];
			break;
		case CSG::VER_2:
			output = tri->corner[2];
			break;
		}
	}
	return output;
}

void ISectTest2(Octree<MPMesh2>* pOctree, bool bInverse)
{
	assert(pOctree);
	std::list<OctreeNode*> leaves;
	GetLeafNodes(pOctree->Root, leaves, NODE_COMPOUND);

	std::set<GS::IndexPair> antiOverlapMap;
	for (auto leaf: leaves)
	{
		unsigned i, j, ni, nj;
		MPMesh2 *mesh1 = pOctree->pMesh[0];
		MPMesh2 *mesh2 = pOctree->pMesh[1];
		MPMesh2::FaceHandle tri1, tri2;
		Vec3d *v0,*v1,*v2, nv,*u0,*u1,*u2,nu,start,end;
		MPMesh::FVIter fvItr;
		int isISect;
		int startT(0), endT(0);
		VertexPos startiT, startjT, endiT, endjT;
		FeitoISectZone::ISVertexItr vP1, vP2;

        int triId[2];
        GS::IndexPair iPair;
        auto &tt1 = leaf->TriangleTable[0];
        auto &tt2 = leaf->TriangleTable[1];
        ni = tt1.size();
        nj = tt2.size();

		for (i = 0; i < ni; i++)
		{
			for (j = 0; j < nj; j++)
			{
				tri1 = tt1[i];
				tri2 = tt2[j];
                triId[0] = tri1.idx();
                triId[1] = tri2.idx();
				GS::MakeIndex(triId, iPair);
				if (antiOverlapMap.find(iPair) != antiOverlapMap.end()) continue;
				else antiOverlapMap.insert(iPair);

				// intersection test main body
				GetCorners(mesh1, tri1, v0, v1, v2);
				GetCorners(mesh2, tri2, u0, u1, u2);

				nv = mesh1->normal(tri1);
				nu = mesh2->normal(tri2);
							
				startT = 0; endT = 0; // return to Zero.

				isISect = TriTriIntersectTest(*v0, *v1, *v2, nv,
					*u0, *u1, *u2, nu, startT, endT, start, end);

				if (isISect < 0) continue;

				FeitoISectZone **si = &mesh1->property(mesh1->SurfacePropHandle, tri1);
				FeitoISectZone **sj = &mesh2->property(mesh2->SurfacePropHandle, tri2);

				if (!*si) *si = new FeitoISectZone(mesh1, mesh2, tri1);
				if (!*sj) *sj = new FeitoISectZone(mesh2, mesh1, tri2);

				if (isISect == 0)
				{
					(*si)->coplanarTris.emplace_back(tri2);
					(*sj)->coplanarTris.emplace_back(tri1);
					continue;
				}

				startiT = VertexPos(startT & 0xffff);
				startjT = VertexPos(startT >> 16);

				endiT = VertexPos(endT & 0xffff);
				endjT = VertexPos(endT >> 16);

				if (IsEqual(start, end))
				{
					// 点相交
					Vec3d point = (start+end)/2;

					// 最后一个参数表示，可能存在两个以上的插入点
                    
                    auto info1 = InsertPoint(*si, VertexPos(startiT), point);
					InsertPoint(*si, VertexPos(endiT), *info1);
					auto info2 = InsertPoint(*sj, VertexPos(startjT), point);
					InsertPoint(*sj, VertexPos(endjT), *info2);

     //               auto info1 = InsertPoint(*si, VertexPos(startiT | endiT), point);
					//auto info2 = InsertPoint(*sj, VertexPos(startjT | endjT), point);
					//if(!info1->vhandle.is_valid() || !info2->vhandle.is_valid())
					//{
					//	(*sj)->coplanarTris.emplace_back(tri1);
					//	return;
					//}

                    auto &id1 = mesh1->property(mesh1->VertexIndexPropHandle, info1->vhandle);
                    auto &id2 = mesh2->property(mesh2->VertexIndexPropHandle, info2->vhandle);

					if (id1 == -1)
                    {
						assert(id2 == -1);
						int id = FeitoISectZone::add_point(point);
						id1 = id;
						id2 = id;
					}
				}
				else
				{
					// 线相交
					double d = OpenMesh::dot(OpenMesh::cross(nv, nu), end-start);
					//if (bInverse) d = -d;
					
					auto info1 = InsertPoint(*si, startiT, start);
					auto info2 = InsertPoint(*si, endiT, end);

                    auto infox1 = InsertPoint(*sj, startjT, start);
					auto infox2 = InsertPoint(*sj, endjT, end);

                    auto &id1 = mesh1->property(mesh1->VertexIndexPropHandle, info1->vhandle);
                    auto &id2 = mesh2->property(mesh2->VertexIndexPropHandle, infox1->vhandle);
					if (id1 == -1)
					{
						assert(id2 == -1);
						int id = FeitoISectZone::add_point(start);
						id1 = id;
						id2 = id;
					}

                    auto &idx1 = mesh1->property(mesh1->VertexIndexPropHandle, info2->vhandle);
                    auto &idx2 = mesh2->property(mesh2->VertexIndexPropHandle, infox2->vhandle);
					if (idx1 == -1)
					{
						assert(idx2 == -1);
						int id = FeitoISectZone::add_point(end);
						idx1 = id;
						idx2 = id;
					}
					
					if (d > 0)
					{
						(*si)->InsertSegment(info1, info2, *sj);
						(*sj)->InsertSegment(infox2, infox1, *si);
					}
					else
					{
						(*si)->InsertSegment(info2, info1, *sj);
						(*sj)->InsertSegment(infox1, infox2, *si);
					}
				}
			}
		}
	}
}

inline int FindMaxIndex(Vec3d& vec)
{
	double a = fabs(vec[0]);
	double b = fabs(vec[1]);
	double c = fabs(vec[2]);

	if (a >= b)
	{
		if (a >= c) return 0;
		else return 2;
	}
	else
	{
		if (b >= c) return 1;
		else return 2;
	}
}

struct Data2D
{
	GEOM_FADE2D::Point2 p2;
	FeitoISectZone::ISVertexItr itr;
};

void FeitoParsingFace1(FeitoISectZone* triangle, Octree<MPMesh2>* pOctree, std::vector<Data2D>& points)
	{
	// 树不为空，存在一个on的节点
	assert(triangle);
	MPMesh2 *pMesh = triangle->thisMesh;

	Vec3d normal = pMesh->normal(triangle->face);
	int mainAxis = FindMaxIndex(normal);
	if (normal[mainAxis] > 0.0)
	{
		triangle->xi = (mainAxis+1)%3;
		triangle->yi = (mainAxis+2)%3;
	}
	else
	{
		triangle->yi = (mainAxis+1)%3;
		triangle->xi = (mainAxis+2)%3;
	}

	unsigned n_vertices = triangle->vertices.size();

	// 3D 转 2D 坐标
	points.reserve(n_vertices+5);
	points.resize(n_vertices);
	size_t count = 0; 
	for (auto vertex = triangle->vertices.begin();
		vertex != triangle->vertices.end(); vertex++, count++)
	{
		vertex->Id_2d = count;
		points[count].p2.set(vertex->coord[triangle->xi], vertex->coord[triangle->yi]);
		points[count].p2.setCustomIndex(count);
		points[count].itr = vertex;
	}
}

inline void CalcLineCoef(FeitoISectZone::SegData& seg, std::vector<Data2D>& infos)
{
	auto &p0 = infos[seg.points[0]->Id_2d].p2;
	auto &p1 = infos[seg.points[1]->Id_2d].p2;
		
	assert(!IsEqual(seg.points[0]->coord, seg.points[1]->coord));

	auto dir = p1-p0;
    OpenMesh::Vec2d n(-dir.y(), dir.x());
    n.normalize();
	double d = n[1]*p0.y()+n[0]*p0.x();
	seg.lineCoef[0] = n[0];
	seg.lineCoef[1] = n[1];
	seg.lineCoef[2] = -d;
}

inline double cross(const OpenMesh::Vec2d& v0, const OpenMesh::Vec2d& v1)
{
	return v0[0]*v1[1]-v0[1]*v1[0];
}

inline bool IsPointInTriangle(const OpenMesh::Vec2d &bc, OpenMesh::Vec2d* v, Relation& rel)
{
	double d = cross(v[1]-v[0], bc-v[0]);
	if (cross(v[2]-v[1], bc-v[1])*d < 0) return false;
	if (cross(v[0]-v[2], bc-v[2])*d < 0) return false;
	return true;
}

bool IsInsideTriangle(FeitoISectZone* triangle, MPMesh2::FaceHandle coTri, const Point2 &bc, Relation &rel)
{
	Vec3d *v0, *v1, *v2;
	auto pMesh = triangle->thatMesh;
	GetCorners(pMesh, coTri, v0, v1, v2);
		
	OpenMesh::Vec2d v[3];
	v[0][0] = (*v0)[triangle->xi]; v[0][1] = (*v0)[triangle->yi];
	v[1][0] = (*v1)[triangle->xi]; v[1][1] = (*v1)[triangle->yi];
	v[2][0] = (*v2)[triangle->xi]; v[2][1] = (*v2)[triangle->yi];

	if (IsPointInTriangle(OpenMesh::Vec2d(bc.x(), bc.y()), v, rel))
	{
		double d = OpenMesh::dot(triangle->thisMesh->normal(triangle->face), pMesh->normal(coTri));
		if (d > 0.0) rel = REL_SAME;
		else rel = REL_OPPOSITE;
		return true;
	}
	else return false;
}

void FeitoParsingFace2(FeitoISectZone* triangle, Octree<MPMesh2>* pOctree, std::vector<Data2D>& points)
{
	BSPSeg tmpSeg;
	std::vector<BSPSeg> segments;
	auto &segs = triangle->segs;
	for (auto &seg: segs)
	{
		CalcLineCoef(seg, points);
		tmpSeg.lineCoef = seg.lineCoef;
		tmpSeg.start = points[seg.points[0]->Id_2d].p2;
		tmpSeg.end = points[seg.points[1]->Id_2d].p2;
		segments.push_back(tmpSeg);
	}
	BSP2D* bsp =  BuildBSP2DNode(segments);

	std::vector<GEOM_FADE2D::Segment2> segList;
	Point2 *p0, *p1;
	for (auto seg = segs.begin(); seg != segs.end(); seg++)
	{
		p0 = &points[seg->points[0]->Id_2d].p2;
		p1 = &points[seg->points[1]->Id_2d].p2;
		segList.emplace_back(*p0, *p1);
	}
	
	GEOM_FADE2D::Fade_2D *dt;
    dt = new GEOM_FADE2D::Fade_2D;
	for (int i = 0; i < points.size(); i++)
		dt->insert(points[i].p2);
	//dt->insert(points[triangle->corner[0]->Id_2d].p2);
	//dt->insert(points[triangle->corner[1]->Id_2d].p2);
	//dt->insert(points[triangle->corner[2]->Id_2d].p2);

	dt->createConstraint(segList, GEOM_FADE2D::CIS_IGNORE_DELAUNAY);
	dt->applyConstraintsAndZones();

	std::vector<GEOM_FADE2D::Triangle2*> vAllTriangles;
	dt->getTrianglePointers(vAllTriangles);

	triangle->thisMesh->delete_face(triangle->face, false);

	auto pMesh = triangle->thisMesh;
	GEOM_FADE2D::Point2 baryCenter2d;
	for (auto triFrag: vAllTriangles)
	{
		baryCenter2d = triFrag->getBarycenter();
		GS::double3 v[3];
		v[0] = Vec3dToDouble3(points[triFrag->getCorner(0)->getCustomIndex()].itr->coord);
		v[1] = Vec3dToDouble3(points[triFrag->getCorner(1)->getCustomIndex()].itr->coord);
		v[2] = Vec3dToDouble3(points[triFrag->getCorner(2)->getCustomIndex()].itr->coord);

		GS::double3x3 mat(GS::double3(1,1,1), v[2]-v[1], v[2]-v[0]);
		if (fabs(GS::determinant(mat)) < 1e-9) continue;

		auto fhandle = pMesh->add_face(
			points[triFrag->getCorner(0)->getCustomIndex()].itr->vhandle,
			points[triFrag->getCorner(1)->getCustomIndex()].itr->vhandle,
			points[triFrag->getCorner(2)->getCustomIndex()].itr->vhandle);

		//assert(fhandle.is_valid());
		if (!fhandle.is_valid()) continue;

		Relation rel = REL_UNKNOWN;
		int id = fhandle.idx();
        auto &cop = triangle->coplanarTris;
		if (cop.size())
        {
			for (auto &coTri: cop)
			{
				if (IsInsideTriangle(triangle, coTri, baryCenter2d, rel))
					break;
			}
        }
		else
		{
			if (bsp)
				rel = BSP2DInOutTest(bsp, &baryCenter2d);
		}
		switch (rel)
		{
		case CSG::REL_INSIDE:
			pMesh->property(pMesh->TopologyInfo, fhandle) = REL_INSIDE;
			break;
		case CSG::REL_OUTSIDE:
			pMesh->property(pMesh->TopologyInfo, fhandle) = REL_OUTSIDE;
			break;
		case CSG::REL_OPPOSITE:
			pMesh->property(pMesh->TopologyInfo, fhandle) = REL_OPPOSITE;
			break;
		case CSG::REL_SAME:
			pMesh->property(pMesh->TopologyInfo, fhandle) = REL_SAME;
			break;
		default:
			break;
		}
	}

	SAFE_RELEASE(dt);
	SAFE_RELEASE(bsp);
}


void Tesselation(Octree<MPMesh2>* pOctree)
{
	std::vector<Data2D> points;
	for (auto &triangle: FeitoISectZone::zone_record)
	{
		points.clear();
		FeitoParsingFace1(triangle, pOctree, points);
		FeitoParsingFace2(triangle, pOctree, points);
	}
}

static int booleanRelation[3][2] = {
	{REL_OUTSIDE | REL_SAME, REL_OUTSIDE},
	{REL_INSIDE | REL_SAME, REL_INSIDE},
	{REL_OUTSIDE | REL_OPPOSITE, REL_INSIDE}
}; 

MPMesh2* Classification(Octree<MPMesh2>* pOctree, MPMesh2** meshes, int op)
{
	MPMesh2* res = new MPMesh2;
	std::queue<MPMesh2::FaceHandle> queue1, queue2;
	std::list<MPMesh2::FaceHandle> faceBuffer;
	for (int main = 0; main < 2; main++)
	{
		auto *pMesh = meshes[main];
		pMesh->garbage_collection();
		MPMesh2::VertexHandle *record = new MPMesh2::VertexHandle[pMesh->n_vertices()];
		for (int i = 0; i < pMesh->n_vertices(); i++) record[i].reset();

		for (auto face_itr = pMesh->faces_sbegin(); 
			face_itr != pMesh->faces_end(); face_itr++)
		{
			if (pMesh->property(pMesh->MarkPropHandle, *face_itr) == 2) continue;
			queue2.push(*face_itr);

			int targetRelation = booleanRelation[op][main];
			int accept, rule;
			while (!queue2.empty())
			{
				if (pMesh->property(pMesh->MarkPropHandle, queue2.front()) == 2)
				{
					queue2.pop();
					continue;
				}
				//assert(pMesh->property(pMesh->TopologyInfo, queue2.front()) != 0);
				queue1.push(queue2.front());
				queue2.pop();
				accept = -1; rule = -1;
				while (!queue1.empty())
				{
					if (pMesh->property(pMesh->MarkPropHandle, queue1.front()) == 2)
					{
						queue1.pop();
						continue;
					}
					auto curFace = queue1.front();
					queue1.pop();
					pMesh->property(pMesh->MarkPropHandle, curFace) = 2;

					faceBuffer.push_back(curFace);

					if (accept == -1)
					{
						int a = pMesh->property(pMesh->TopologyInfo, curFace);
						if (a)
						{
							if (a & targetRelation) accept = 1;
							else accept = 0;
							rule = a;
						}
					}

					MPMesh2::FaceFaceIter ffItr = pMesh->ff_iter(curFace);
					int *markPtr;
					for (int i = 0; i < 3; i++, ffItr++)
					{
						if (!ffItr.is_valid())
							break;

						markPtr = &(pMesh->property(pMesh->MarkPropHandle, *ffItr));
						if (*markPtr != 2)
						{
							int b = pMesh->property(pMesh->TopologyInfo, *ffItr);
							if (rule > -1 && b)
							{
								if (rule == b && *markPtr != 1)
								{
									*markPtr = 1; // queued
									queue1.push(*ffItr);
								}
								else if (rule != b && *markPtr != 4) // seed
								{
									*markPtr = 4; // seed
									queue2.push(*ffItr);
								}
							}
							else
							{
								if (*markPtr != 1)
								{
									*markPtr = 1; // queued
									queue1.push(*ffItr);
								}
							}
						}
					}
				} // queue1

				assert(accept != -1);
				//static int count = 0;
				//int count2 = count;
				if (accept == 1)
				{
					for (auto itr:faceBuffer)
					{
						//if (!(count2--)) break;;
						//debug
						//auto dn = pMesh->n_faces();
						auto fvItr = pMesh->fv_begin(itr);

						MPMesh2::VertexHandle *vhandle3[3];
						for (int i = 0; i < 3; i++,fvItr++)
						{
							// 检查是否已被添加
							int global_id = pMesh->property(pMesh->VertexIndexPropHandle, *fvItr).value;
							if (global_id != -1)
							{
								vhandle3[i] = &(FeitoISectZone::point(global_id).used);
								if (!vhandle3[i]->is_valid())
									*vhandle3[i] = res->add_vertex(pMesh->point(*fvItr));
							
								record[fvItr->idx()] = *(vhandle3[i]);
							}
							else
							{
								vhandle3[i] = &record[fvItr->idx()];
								if (!record[fvItr->idx()].is_valid())
									*vhandle3[i] = res->add_vertex(pMesh->point(*fvItr));
							}
						}

						MPMesh2::FaceHandle fhandle;
						if (main*op == 2)
							fhandle = res->add_face(*vhandle3[2], *vhandle3[1], *vhandle3[0]);
						else 
							fhandle = res->add_face(*vhandle3[0], *vhandle3[1], *vhandle3[2]);
						//assert(fhandle.is_valid());
					}
				}
				faceBuffer.clear();
			//count++;
			//break;
			} // queue2
		}
		delete [] record;
	}
	return res;
}


// 0 Union:
// 1 Intersect
// 2 Diff 

extern "C" CSG_API MPMesh2* BooleanOperationFeito(MPMesh2* mesh1, MPMesh2* mesh2, int op, HANDLE hd)
{
    _output = hd;
    StdOutput("Hi");

    MPMesh2* res = nullptr;
    MPMesh2 *arrMesh[2] = {mesh1, mesh2};
    Octree<MPMesh2>* pOctree = BuildOctree2(arrMesh, 2);
    ISectTest2(pOctree, (op==2)?true:false);
    Tesselation(pOctree);
    res = Classification(pOctree, arrMesh, op);

    FeitoISectZone::clear();
    delete pOctree;

	//static int count = 0;
	//count ++;
	//if (count == 2)
	//{
	////res = arrMesh[1];
	//}
	//res->garbage_collection();
	size_t n = res->n_vertices();
	res->verticesList = new Vec3d[n];

	for (auto vItr = res->vertices_begin();
		vItr != res->vertices_end(); vItr ++)
		res->verticesList[vItr->idx()] = res->point(*vItr);

	res->update_normals();
	//GS::BaseMesh *p = Convert2BaseMesh2(res);
	//delete res;
	//res = Convert2MPMesh2(p);
	//delete p;

    return res;
}

inline OpenMesh::Vec3d convert_double3(GS::double3& vec)
{
	return OpenMesh::Vec3d(vec.x, vec.y, vec.z);
}

extern "C" CSG_API MPMesh2* Convert2MPMesh2(const GS::BaseMesh* pMesh)
{
    if (!pMesh) return nullptr;

	MPMesh2 *res = new MPMesh2(pMesh);

    int n = (int)pMesh->VertexCount();
    res->BBox.Clear();
	MPMesh::VertexHandle vhandle;
	res->verticesList = new Vec3d[n];
    for (int i = 0; i < n; i++)
	{
		res->verticesList[i] = convert_double3(pMesh->Vertex(i).pos);
        vhandle = res->add_vertex(res->verticesList[i]);
		//res->property(res->PointInOutTestPropHandle, vhandle) = 0;
		res->BBox.IncludePoint(res->point(vhandle));
	}

	std::vector<MPMesh::VertexHandle>  face_vhandles(3);
    n = (int)pMesh->PrimitiveCount();
    for (int i = 0; i < n; i++)
    {
        auto &info = pMesh->TriangleInfo(i);
		face_vhandles[0] = res->vertex_handle(info.VertexId[0]);
		face_vhandles[1] = res->vertex_handle(info.VertexId[1]);
		face_vhandles[2] = res->vertex_handle(info.VertexId[2]);
        res->add_face(face_vhandles);
    }

	res->update_normals();

	return res;
}

extern "C" CSG_API GS::BaseMesh* Convert2BaseMesh2(MPMesh2* mesh)
{
    if (!mesh) return nullptr;
    GS::BaseMesh* res = new GS::BaseMesh;

    for (auto face = mesh->faces_begin(); face != mesh->faces_end(); face++)
    {
		GS::double3 v[3];
		Vec3d *v0, *v1, *v2;
		GetCorners(mesh, static_cast<MPMesh::FaceHandle>(*face), v0, v1, v2);
		v[0] = Vec3dToDouble3(*v0);
		v[1] = Vec3dToDouble3(*v1);
		v[2] = Vec3dToDouble3(*v2);
		res->AddTriangle(v);
    }

    return res;
}


}
