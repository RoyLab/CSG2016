#include "precompile.h"
#include <CGAL\Plane_3.h>
#include <CGAL\intersection_3.h> 
#include <algorithm>

#include <xmath.h>
#include <xlogger.h>
#include "csgdefs.h"
//#include "mathext.h"
#include <macros.h>
//#include "MyAlgorithm.h"
#include "Octree.h"
//#include "MyMesh.h"
//#include "AlgUserData.h"

namespace Boolean
{
    typedef Depeck K;
    typedef K::Vector_3 Vec3d;
    typedef K::Point_3  Point3d;
    typedef K::Ray_3    Ray;

    enum RayTriRelation
    {
        RTR_NONE,
        RTR_INTERSECT,
        RTR_ERROR,
        RTR_ON
    };

    static int _ray_count = 0;

    template <class FT>
    static int FindFirstNode(FT tx0, FT ty0, FT tz0, FT txm, FT tym, FT tzm)
    {
        unsigned char answer = 0; // initialize to 00000000
                                  // select the entry plane and set bits
        if (tx0 > ty0) {
            if (tx0 > tz0) { // PLANE YZ
                if (tym < tx0) answer |= 2; // set bit at position 1
                if (tzm < tx0) answer |= 1; // set bit at position 0 			
                return (int)answer;
            }
        }
        else if (ty0 > tz0) { // PLANE XZ
            if (txm < ty0) answer |= 4; // set bit at position 2
            if (tzm < ty0) answer |= 1; // set bit at position 0
            return (int)answer;

        }
        // PLANE XY
        if (txm < tz0) answer |= 4;   // set bit at position 2
        if (tym < tz0) answer |= 2;   // set bit at position 1
        return (int)answer;
    }

    template <class FT>
    static int GetNextNode(FT txm, int x, FT tym, int y, FT tzm, int z)
    {
        if (txm < tym) {
            if (txm < tzm) { return x; }  // YZ plane
        }
        else {
            if (tym < tzm) { return y; } // XZ plane
        }
        return z; // XY plane;
    }

    struct RayCastInfo
    {
        Ray ray;
        int rayId, meshId;
        int nCross = 0;
    };

    static RayTriRelation RayFaceTest(K::Ray_3 &ray, K::Triangle_3& triangle)
    {
        if (triangle.is_degenerate())
            return RTR_NONE;

        bool result = CGAL::do_intersect(ray, triangle);
        if (!result) return RTR_NONE;

        if (triangle.has_on(ray.source()))
            return RTR_ON;

        for (size_t i = 0; i < 3; i++)
            if (CGAL::coplanar(ray.source(), ray.point(1), triangle.vertex(i), triangle.vertex((i + 1) % 3)))
                return RTR_ERROR;

        return RTR_INTERSECT;
    }

    static RayTriRelation RayFaceTest(K::Ray_3 &ray, Triangle* triangle)
    {
        K::Triangle_3 cgal_tri = convertToCGALTriangle<Depeck>(triangle);
        return RayFaceTest(ray, cgal_tri);
    }


    static void RayCastThroughNode(Octree::Node* pNode, RayCastInfo& rayInfo)
    {
        if (pNode == NULL) return;

        auto trianglesItr = pNode->triTable.find(rayInfo.meshId);
        if (trianglesItr == pNode->triTable.end()) return;

        auto &triangles = trianglesItr->second;
        unsigned n = triangles.size();

        for (unsigned i = 0; i < n; i++)
        {
            auto fhandle = (triangles)[i];
            if (fhandle->octree_mark == rayInfo.rayId) continue;

            fhandle->octree_mark = rayInfo.rayId;
            RayTriRelation ret = RayFaceTest(rayInfo.ray, fhandle);

            switch (ret)
            {
            case RTR_INTERSECT:
                rayInfo.nCross += 1;
                break;
            case RTR_ON:
            case RTR_ERROR:
                throw ret;
            default:
                break;
            }
        }
    }

    template <class FT>
    static void ProcessSubNode(FT tx0, FT ty0, FT tz0, FT tx1, FT ty1, FT tz1,
        Octree::Node* pNode, int& a, RayCastInfo* rayInfo)
    {
        if (pNode == nullptr) return;

        FT txm, tym, tzm;
        int currNode;
        if (tx1 < 0 || ty1 < 0 || tz1 < 0)
            return;
        if (pNode->type != NODE_MIDSIDE)
        {
            RayCastThroughNode(pNode, *rayInfo);
            return;
        }

        txm = 0.5*(tx0 + tx1);
        tym = 0.5*(ty0 + ty1);
        tzm = 0.5*(tz0 + tz1);
        currNode = FindFirstNode(tx0, ty0, tz0, txm, tym, tzm);
        do {
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

    static void RayTraverse(Octree* pOctree, RayCastInfo &rayInfo, const cyPointT& source)
    {
        assert(pOctree);
        const XR::BoundingBox &RootBBox = pOctree->getRoot()->bbox;
        cyPointT bboxSize = RootBBox.diagonal<cyPointT>();
        cyPointT source_copy = source;
        cyPointT dir = cyPointT(
            randf() + 0.001,
            randf() + 0.001,
            randf() + 0.001
        ); // no zero

        int a = 0;
        for (int i = 0; i < 3; ++i)
        {
            if (dir[i] < 0.0)
            {
                source_copy[i] = bboxSize[i] - source[i];
                dir[i] = -dir[i];
                a |= (1 << (2 - i));
            }
        }

        rayInfo.ray = Ray(
            Point3d(source_copy[0], source_copy[1], source_copy[2]),
            K::Direction_3(dir[0], dir[1], dir[2])
        );

        cyPointT t0 = (RootBBox.min<cyPointT>() - source_copy);
        cyPointT t1 = (RootBBox.max<cyPointT>() - source_copy);

        //Point3d cgal_source(source_copy[0], source_copy[1], source_copy[2]);
        //Vec3d t0 = (RootBBox.min<Point3d>() - cgal_source);
        //Vec3d t1 = (RootBBox.max<Point3d>() - cgal_source);

        if (std::max(t0[0], std::max(t0[1], t0[2])) < std::min(t1[0], std::min(t1[1], t1[2])))
        {
            ProcessSubNode(
                t0[0] / dir[0], t0[1] / dir[1], t0[2] / dir[2],
                t1[0] / dir[0], t1[1] / dir[1], t1[2] / dir[2],
                pOctree->getRoot(), a, &rayInfo);
        }
        else
        {
            //XLOG_DEBUG << "No test.";
        }
    }

    static Ray randRay(Point3d& pt)
    {
        Vec3d dir(randf(), randf(), randf());
        return Ray(pt, dir);
    }

    Relation PolyhedralInclusionTest(cyPointT& point, Octree* pOctree,
        std::vector<RegularMesh*>& pMesh, unsigned meshId, bool IsInverse, int *cross)
    {
        if (pMesh[meshId]->bbox().has_on_unbounded_side(point))
        {
            if (!IsInverse) return REL_OUTSIDE;
            else return REL_INSIDE;
        }

        auto bbox = pOctree->getRoot()->bbox;
        static int rayId = 1;

        RayCastInfo rayInfo;
        rayInfo.meshId = meshId;

        bool isValid = false;
        int n_try = 0;
        while (!isValid)
        {
            rayInfo.rayId = rayId++;
            rayInfo.nCross = 0;
            isValid = true;
            try
            {
                RayTraverse(pOctree, rayInfo, point);
            }
            catch (RayTriRelation result)
            {
                switch (result)
                {
                case RTR_ERROR:
                    XLOG_DEBUG << "Regen direction " << ++n_try;
                    isValid = false;
                    break;
                case RTR_ON:
                    return REL_ON_BOUNDARY;
                default:
                    ReportError("");
                    return REL_NOT_AVAILABLE;
                    break;
                }
            }
        }

        if (cross)
        {
            *cross = rayInfo.nCross;
        }

        if (IsInverse)
        {
            if (rayInfo.nCross % 2 == 1)
                return REL_OUTSIDE;
            else return REL_INSIDE;
        }
        else
        {
            if (rayInfo.nCross % 2 == 1)
                return REL_INSIDE;
            else return REL_OUTSIDE;
        }
    }

}