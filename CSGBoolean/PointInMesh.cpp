#include <CGAL\Plane_3.h>
#include <CGAL\intersection_3.h>
#include "csgdefs.h"

#include "MyAlgorithm.h"
#include "Octree.h"
#include "MyMesh.h"
#include "AlgUserData.h"

namespace CSG
{
    typedef K::Vector_3 Vec3d;
    typedef K::Point_3  Point3d;

    enum RayTriRelation
    {
        RTR_NONE,
        RTR_INTERSECT,
        RTR_ERROR,
        RTR_ON
    };

    static int _ray_count = 0;

    static int FindFirstNode(double tx0, double ty0, double tz0, double txm, double tym, double tzm)
    {
        unsigned char answer = 0; // initialize to 00000000
        // select the entry plane and set bits
        if (tx0 > ty0){
            if (tx0 > tz0){ // PLANE YZ
                if (tym < tx0) answer |= 2; // set bit at position 1
                if (tzm < tx0) answer |= 1; // set bit at position 0 			
                return (int)answer;
            }
        }
        else if (ty0 > tz0){ // PLANE XZ
            if (txm < ty0) answer |= 4; // set bit at position 2
            if (tzm < ty0) answer |= 1; // set bit at position 0
            return (int)answer;

        }
        // PLANE XY
        if (txm < tz0) answer |= 4;   // set bit at position 2
        if (tym < tz0) answer |= 2;   // set bit at position 1
        return (int)answer;
    }

    static int GetNextNode(double txm, int x, double tym, int y, double tzm, int z)
    {
        if (txm < tym){
            if (txm < tzm){ return x; }  // YZ plane
        }
        else{
            if (tym < tzm){ return y; } // XZ plane
        }
        return z; // XY plane;
    }

    struct RayCastInfo
    {
        K::Ray_3 ray;
        int rayId, meshId;
        int nCross = 0;
        K::Plane_3 splane;
    };

    //static inline bool PointWithPlane(const Vec3d& normal, double distance, const Vec3d& pos)
    //{
    //    double dist = pos * normal + distance;
    //    if (dist > 0)
    //        return true;
    //    return false;
    //}

    //inline bool PointWithPlane(const Vec3d& normal, const Vec3d& v0, const Vec3d& pos)
    //{
    //    double dist = normal * (pos - v0);
    //    if (dist > 0)
    //        return true;
    //    return false;
    //}

    //static inline bool SignDeterminant(const Point3d& v0, const Point3d& v1, const Point3d& v2, const Point3d& v3)
    //{
    //    return CGAL::orientation(v0, v1, v2, v3);
    //}

    //static inline void normalize(Vec3d& vec)
    //{
    //    vec = vec / sqrt(vec.squared_length());
    //}

    static RayTriRelation RayFaceTest(K::Ray_3 &ray, K::Triangle_3& triangle)
    {
        if (triangle.has_on(ray.source()))
            return RTR_ON;

        bool result = CGAL::do_intersect(ray, triangle);
        if (!result) return RTR_NONE;

        for (size_t i = 0; i < 3; i++)
            if (CGAL::coplanar(ray.source(), ray.point(1), triangle.vertex(i), triangle.vertex((i + 1) % 3)))
                return RTR_ERROR;

        return RTR_INTERSECT;
    }

    static void RayCastThroughNode(Octree::Node* pNode, RayCastInfo* rayInfo)
    {
        if (pNode == NULL) return;

        auto trianglesItr = pNode->triTable.find(rayInfo->meshId);
        if (trianglesItr == pNode->triTable.end()) return;

        auto &triangles = trianglesItr->second;
        unsigned n = triangles->size();

        for (unsigned i = 0; i < n; i++)
        {
            auto fhandle = (*triangles)[i];
            if (fhandle->mark == rayInfo->rayId) continue;
            
            fhandle->mark = rayInfo->rayId;
            RayTriRelation ret = RayFaceTest(rayInfo->ray, fhandle->data->triangle);

            switch (ret)
            {
            case RTR_INTERSECT:
                rayInfo->nCross ^= 1;
                break;
            case RTR_ON:
            case RTR_ERROR:
                throw ret;
            default:
                break;
            }
        }
    }

    static RayTriRelation ProcessSubNode(double tx0, double ty0, double tz0, double tx1, double ty1, double tz1,
        Octree::Node* pNode, int& a, RayCastInfo* rayInfo)
    {
        if (pNode == nullptr) return;

        double txm, tym, tzm;
        int currNode;
        if (tx1 < 0 || ty1 < 0 || tz1 < 0)
            return;
        if (pNode->type != NODE_MIDSIDE)
        {
            RayCastThroughNode(pNode, rayInfo);
            return;
        }

        txm = 0.5*(tx0 + tx1);
        tym = 0.5*(ty0 + ty1);
        tzm = 0.5*(tz0 + tz1);
        currNode = FindFirstNode(tx0, ty0, tz0, txm, tym, tzm);
        do{
            if (rayInfo->nCross < 0) return;

            switch (currNode)
            {
            case 0: {
                ProcessSubNode(tx0, ty0, tz0, txm, tym, tzm, &pNode->pChildren[a], a, rayInfo);
                currNode = GetNextNode(txm, 4, tym, 2, tzm, 1);
                break; }
            case 1: {
                ProcessSubNode(tx0, ty0, tzm, txm, tym, tz1, &pNode->pChildren[1 ^ a], a, rayInfo);
                currNode = GetNextNode(txm, 5, tym, 3, tz1, 8);
                break; }
            case 2: {
                ProcessSubNode(tx0, tym, tz0, txm, ty1, tzm, &pNode->pChildren[2 ^ a], a, rayInfo);
                currNode = GetNextNode(txm, 6, ty1, 8, tzm, 3);
                break; }
            case 3: {
                ProcessSubNode(tx0, tym, tzm, txm, ty1, tz1, &pNode->pChildren[3 ^ a], a, rayInfo);
                currNode = GetNextNode(txm, 7, ty1, 8, tz1, 8);
                break; }
            case 4: {
                ProcessSubNode(txm, ty0, tz0, tx1, tym, tzm, &pNode->pChildren[4 ^ a], a, rayInfo);
                currNode = GetNextNode(tx1, 8, tym, 6, tzm, 5);
                break; }
            case 5: {
                ProcessSubNode(txm, ty0, tzm, tx1, tym, tz1, &pNode->pChildren[5 ^ a], a, rayInfo);
                currNode = GetNextNode(tx1, 8, tym, 7, tz1, 8);
                break; }
            case 6: {
                ProcessSubNode(txm, tym, tz0, tx1, ty1, tzm, &pNode->pChildren[6 ^ a], a, rayInfo);
                currNode = GetNextNode(tx1, 8, ty1, 8, tzm, 7);
                break; }
            case 7: {
                ProcessSubNode(txm, tym, tzm, tx1, ty1, tz1, &pNode->pChildren[7 ^ a], a, rayInfo);
                currNode = 8;
                break; }
            }
        } while (currNode < 8);

    }

    static void RayTraverse(Octree* pOctree, RayCastInfo* rayInfo)
    {
        assert(pOctree);
        AABBmp &RootBBox = pOctree->BoundingBox;
        static unsigned count = 0; // for every new raytraverse a new Id
        count++;

        int a = 0;
        if (rayInfo->dir[0] < 0.0)
        {
            rayInfo->pt[0] = RootBBox.Size()[0] - rayInfo->pt[0];
            rayInfo->dir[0] = -rayInfo->dir[0];
            a |= 4;
        }
        if (rayInfo->dir[1] < 0.0)
        {
            rayInfo->pt[1] = RootBBox.Size()[1] - rayInfo->pt[1];
            rayInfo->dir[1] = -rayInfo->dir[1];
            a |= 2;
        }
        if (rayInfo->dir[2] < 0.0)
        {
            rayInfo->pt[2] = RootBBox.Size()[2] - rayInfo->pt[2];
            rayInfo->dir[2] = -rayInfo->dir[2];
            a |= 1;
        }
        Vec3d t0 = (RootBBox.Min() - rayInfo->pt) / rayInfo->dir;
        Vec3d t1 = (RootBBox.Max() - rayInfo->pt) / rayInfo->dir;
        if (std::max(t0[0], GS::max(t0[1], t0[2])) < GS::min(t1[0], GS::min(t1[1], t1[2])))
            ProcessSubNode(t0[0], t0[1], t0[2], t1[0], t1[1], t1[2], pOctree->Root, a, pMesh, rayInfo, count);
    }

    Relation PolyhedralInclusionTest(Vec3d& point, Octree* pOctree, unsigned meshId, bool IsInverse)
    {
        if (!IsInverse && !pOctree->pMesh[meshId]->BBox.IsInBox(point))
            return REL_OUTSIDE;
        else if (IsInverse && !pOctree->pMesh[meshId]->BBox.IsInBox(point))
            return REL_INSIDE;

        AABBmp &bbox = pOctree->Root->BoundingBox;

        RayCastInfo rayInfo;
        rayInfo.nCross = 0;
        rayInfo.pt = point;
        rayInfo.et = bbox.Max() + Vec3d(1, 1, 1);
        rayInfo.dir = rayInfo.et - point;
        normalize(rayInfo.dir);

        Vec3d edge = rayInfo.et - Vec3d(0.0, bbox.Size()[1] * 0.5, 0.0) - point;
        Vec3d norm = cross(rayInfo.dir, edge);
        normalize(norm);
        GS::double3 tmp1 = Vec3dToDouble3(norm);
        GS::double3 tmp2 = Vec3dToDouble3(point);
        rayInfo.splane = GS::Plane<double>(tmp1, tmp2);
        RayTraverse(pOctree, pOctree->pMesh[meshId], &rayInfo);

        if (IsInverse)
        {
            switch (rayInfo.nCross)
            {
            case -1:    assert(0); return REL_OPPOSITE;
            case -2:    assert(0); return REL_SAME;
            case 1:     return REL_OUTSIDE;
            default:    return REL_INSIDE;
            }
        }

        switch (rayInfo.nCross)
        {
        case -1:    return REL_SAME;
        case -2:    return REL_OPPOSITE;
        case 1:     return REL_INSIDE;
        default:    return REL_OUTSIDE;
        }
    }

    Relation pointInPolyhedron(CGAL::Point_3<K>& p, MyMesh* mesh, Octree* pOctree)
    {

    }

}