#include <CGAL\Plane_3.h>
#include <CGAL\intersection_3.h>
#include <algorithm>

#include "csgdefs.h"
#include "mathext.h"
#include "macroutil.h"
#include "MyAlgorithm.h"
#include "Octree.h"
#include "MyMesh.h"
#include "AlgUserData.h"

namespace CSG
{
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
        Ray ray;
        int rayId, meshId;
        int nCross = 0;
    };

    static RayTriRelation RayFaceTest(K::Ray_3 &ray, K::Triangle_3& triangle)
    {
        bool result = CGAL::do_intersect(ray, triangle);
        if (!result) return RTR_NONE;

        if (triangle.has_on(ray.source()))
            return RTR_ON;

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

    static void ProcessSubNode(double tx0, double ty0, double tz0, double tx1, double ty1, double tz1,
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
        Cube_3 &RootBBox = pOctree->getRoot()->bbox;
        Vec3d  bboxSize = RootBBox.max() - RootBBox.min();

        double dir[3] = { rayInfo->ray.direction().dx(), rayInfo->ray.direction().dy(), rayInfo->ray.direction().dz() };
        double source[3] = { rayInfo->ray.source()[0], rayInfo->ray.source()[1], rayInfo->ray.source()[2] };

        int a = 0;
        if (dir[0] < 0.0)
        {
            source[0] = bboxSize[0] - rayInfo->ray.source()[0];
            dir[0] = -rayInfo->ray.direction().dx();
            a |= 4;
        }
        if (dir[1] < 0.0)
        {
            source[1] = bboxSize[1] - rayInfo->ray.source()[1];
            dir[1] = -rayInfo->ray.direction().dy();
            a |= 2;
        }
        if (dir[2] < 0.0)
        {
            source[2] = bboxSize[2] - rayInfo->ray.source()[2];
            dir[2] = -rayInfo->ray.direction().dz();
            a |= 1;
        }
        rayInfo->ray = Ray(Point3d(source[0], source[1], source[2]), Vec3d(dir[0], dir[1], dir[2]));

        vec3 t0 = { (RootBBox.xmin() - rayInfo->ray.source()[0]) / rayInfo->ray.direction().dx(),
            (RootBBox.ymin() - rayInfo->ray.source()[1]) / rayInfo->ray.direction().dy(),
            (RootBBox.zmin() - rayInfo->ray.source()[2]) / rayInfo->ray.direction().dz() };

        vec3 t1 = { (RootBBox.xmax() - rayInfo->ray.source()[0]) / rayInfo->ray.direction().dx(),
            (RootBBox.ymax() - rayInfo->ray.source()[1]) / rayInfo->ray.direction().dy(),
            (RootBBox.zmax() - rayInfo->ray.source()[2]) / rayInfo->ray.direction().dz() };

        if (std::max(t0[0], std::max(t0[1], t0[2])) < std::min(t1[0], std::min(t1[1], t1[2])))
            ProcessSubNode(t0[0], t0[1], t0[2], t1[0], t1[1], t1[2], pOctree->getRoot(), a, rayInfo);
    }

    static Ray randRay(Point3d& pt)
    {
        Vec3d dir(randf(), randf(), randf());
        return Ray(pt, dir);
    }

    Relation PolyhedralInclusionTest(Point3d& point, Octree* pOctree, std::vector<MyMesh*>& pMesh, unsigned meshId, bool IsInverse)
    {
        if (pMesh[meshId]->get_bbox_cube().has_on_unbounded_side(point))
        {
            if (!IsInverse) return REL_OUTSIDE;
            else return REL_INSIDE;
        }

        auto bbox = pOctree->getRoot()->bbox;
        static int rayId = MARK_BEGIN + 1;

        RayCastInfo rayInfo;
        rayInfo.nCross = 0;
        rayInfo.meshId = meshId;

        bool isValid = false;
        while (!isValid)
        {
            rayInfo.rayId = rayId++;
            rayInfo.ray = randRay(point);
            isValid = true;
            try
            {
                RayTraverse(pOctree, &rayInfo);
            }
            catch (RayTriRelation result)
            {
                switch (result)
                {
                case RTR_ERROR:
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

        if (IsInverse)
        {
            if (rayInfo.nCross == 1)
                return REL_OUTSIDE;
            else return REL_INSIDE;
        }
        else
        {
            if (rayInfo.nCross == 1)
                return REL_INSIDE;
            else return REL_OUTSIDE;
        }
    }

}