#include "csgdefs.h"
#include "ItstAlg.h"
#include "plane_reps.h"

namespace CSG
{

    ItstAlg::ItstAlg(std::vector<MyMesh*>* meshes)
    {
        mp_meshList = meshes;
    }


    ItstAlg::~ItstAlg()
    {
        SAFE_DELETE(mp_adjGraph);
    }


    void ItstAlg::doIntersection(std::vector<Octree::Node*>& intersectLeaves)
    {


        MeshIdTriIdMap antiOverlapMap;
        antiOverlapMap.max_load_factor(0.6);
        mp_adjGraph = new AdjacentGraph(mp_meshList->size());
        auto &meshList = *mp_meshList;

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
                    MyMesh* meshes[2] = { meshList[meshId[0]], meshList[meshId[1]] };

                    // 这里假设map的遍历是保序的，因此meshIdPair自动的分为大小
                    IndexPair meshIdPair;
                    MakeIndex(meshId, meshIdPair);

                    TriIdSet* antiOverlapSet = nullptr;
                    auto searchRes = antiOverlapMap.find(meshIdPair);
                    if (searchRes == antiOverlapMap.end())
                    {
                        antiOverlapSet = new TriIdSet;
                        antiOverlapSet->max_load_factor(0.6);
                        antiOverlapMap.emplace(meshIdPair, antiOverlapSet);
                    }
                    else antiOverlapSet = searchRes->second;

                    for (MyMesh::Face_handle fh0 : *meshItr[0]->second)
                    {
                        for (MyMesh::Face_handle fh1 : *meshItr[1]->second)
                        {
                            uint32_t triId[2] = { fh0->id(), fh1->id() };
                            IndexPair triIdPair;
                            MakeIndex(triId, triIdPair);

                            if (antiOverlapSet->find(triIdPair) != antiOverlapSet->end())
                                continue;

                            antiOverlapSet->insert(triIdPair);

                            if (IntersectionTest(fh0, fh1, antiOverlapSet))
                                mp_adjGraph->setValue(meshId[0], meshId[1], true);
                        }
                    }
                }
            }
        }
    }

    enum Sign
    {
        UNKOWN = 0,
        INTERSECT_ON_LINE,
        INTERSECT_ON_POINT,
        COPLANAR,
        NOT_INTERSECT
    };

    enum PosTag
    {
        NONE = -1, INNER = 0x00,
        EDGE_0 = 0x01, EDGE_1 = 0x02, EDGE_2 = 0x04,
        VER_0 = 0x08, VER_1 = 0x10, VER_2 = 0x20
    };


    template <class _R>
    struct TriTriIsectResult
    {
        PosTag taga[2], tagb[2];

    };

    //template <class _R>
#include "csgdefs.h"
    typedef K _R;

    Sign tri_tri_intersect(const CGAL::Triangle_3<_R>& a, const CGAL::Vector_3<_R>& na,
        const CGAL::Triangle_3<_R>& b, const CGAL::Vector_3<_R>& nb, TriTriIsectResult<_R>* result = nullptr)
    {
        using namespace CGAL;
        typedef _R::FT FT;

        FT d1 = -dot(a.vertex(0), na);
        FT db[3];
        for (size_t i = 0; i < 3; i++)
            db[0] = dot(na, b.vertex(i)) + d1;

        if ((db[0] == db[1]) && (db[1] == db[2]))
        {
            if (db[0] == 0) return COPLANAR;
            else return NOT_INTERSECT;
        }

        FT d2 = -dot(b.vertex(0), nb);
        int a1, b1, c;
        FT da[3];
        for (size_t i = 0; i < 3; i++)
            da[0] = dot(nb, a.vertex(i)) + d2;

        if ((da[0] == da[1]) && (da[1] == da[2]))
        {
            if (da[0] == 0) return COPLANAR;
            else return NOT_INTERSECT;
        }

        /* compute direction of intersection line */
        Vector_3<_R> LineDir = CGAL::cross_product(na, nb);

        /* compute and index to the largest component of D */
        int maxIndex = max_coord(LineDir);


        //int type00(0), type01(0), type10(0), type11(0);

        ///* compute interval for triangle 1 */
        //double isect1[2];
        //Vec3d  isectpointA1, isectpointA2;
        //int res1 = compute_intervals_isectline(v0, v1, v2, index, dv, sdv, &isect1[0], &isect1[1], isectpointA1, isectpointA2, type00, type01);

        ///* compute interval for triangle 2 */
        //double isect2[2];
        //Vec3d  isectpointB1, isectpointB2;
        //int res2 = compute_intervals_isectline(u0, u1, u2, index, du, sdu, &isect2[0], &isect2[1], isectpointB1, isectpointB2, type10, type11);

        //res2 <<= 16;
        //type10 <<= 16;
        //type11 <<= 16;

        //int smallest1 = sort2(isect1[0], isect1[1]);
        //int smallest2 = sort2(isect2[0], isect2[1]);
        //if (isect1[1]<isect2[0] || isect2[1]<isect1[0]) // 要不要用epsf?
        //    return -1;

        //    /* at this point, we know that the triangles intersect */
        //    double dstart = isect1[0] - isect2[0];
        //    if (dstart > EPSF)
        //    {
        //        if (smallest1 == 0)
        //        {
        //            start = isectpointA1;
        //            startType ^= type00;
        //        }
        //        else
        //        {
        //            start = isectpointA2;
        //            startType ^= type01;
        //        }

        //        if (res2) startType ^= res2;
        //    }
        //    else if (dstart < -EPSF)
        //    {
        //        if (smallest2 == 0)
        //        {
        //            start = isectpointB1;
        //            startType ^= type10;
        //        }
        //        else
        //        {
        //            start = isectpointB2;
        //            startType ^= type11;
        //        }
        //        if (res1) startType ^= res1; // 这里可能会出问题
        //    }
        //    else
        //    {
        //        if (smallest1 == 0)
        //        {
        //            start = isectpointA1;
        //            startType ^= type00;
        //        }
        //        else
        //        {
        //            start = isectpointA2;
        //            startType ^= type01;
        //        }

        //        if (smallest2 == 0)
        //        {
        //            start = (isectpointB1 + start) / 2.0;
        //            startType ^= type10;
        //        }
        //        else
        //        {
        //            start = (isectpointB2 + start) / 2.0;
        //            startType ^= type11;
        //        }
        //    }

        //    double dend = isect2[1] - isect1[1];
        //    if (dend > EPSF)
        //    {
        //        if (smallest1 == 0)
        //        {
        //            end = isectpointA2;
        //            endType ^= type01;
        //        }
        //        else
        //        {
        //            end = isectpointA1;
        //            endType ^= type00;
        //        }
        //        if (res2) endType ^= res2;
        //    }
        //    else if (dend < -EPSF)
        //    {
        //        if (smallest2 == 0)
        //        {
        //            end = isectpointB2;
        //            endType ^= type11;
        //        }
        //        else
        //        {
        //            end = isectpointB1;
        //            endType ^= type10;
        //        }
        //        if (res1) endType ^= res1;
        //    }
        //    else
        //    {
        //        if (smallest1 == 0)
        //        {
        //            end = isectpointA2;
        //            endType ^= type01;
        //        }
        //        else
        //        {
        //            end = isectpointA1;
        //            endType ^= type00;
        //        }

        //        if (smallest2 == 0)
        //        {
        //            end = (isectpointB2 + end) / 2.0;
        //            endType ^= type11;
        //        }
        //        else
        //        {
        //            end = (isectpointB1 + end) / 2.0;
        //            endType ^= type10;
        //        }
        //    }

        return 1;
    }


    void checkManifoldEdge()
    {

    }


    bool ItstAlg::IntersectionTest(FH fh0, FH fh1, TriIdSet* overlaps)
    {
        K::Triangle_3& t0 = fh0->data->triangle;
        K::Triangle_3& t1 = fh1->data->triangle;

        // 这里不需要归一化
        K::Vector_3 n0 = CGAL::normal(t0.vertex(0), t0.vertex(1), t0.vertex(2));
        K::Vector_3 n1 = CGAL::normal(t0.vertex(0), t0.vertex(1), t0.vertex(2));

        assert(std::abs(n0.squared_length() - 1) > 1e-6);
        assert(std::abs(n1.squared_length() - 1) > 1e-6);

        TriTriIsectResult<K> result;
        Sign sign = tri_tri_intersect(t0, n0, t1, n1, &result);

        if (sign == NOT_INTERSECT || sign == INTERSECT_ON_POINT || sign == COPLANAR)
            return false;

        ItstTriangle*& it0 = fh0->data->itstTri;
        ItstTriangle*& it1 = fh1->data->itstTri;

        if (!it0) it0 = new ItstTriangle;
        if (!it1) it1 = new ItstTriangle;

        checkManifoldEdge(fh0, fh1, &result, antiOverlapSet);
        setupPonits(fh0, fh1, result);
    }


    //int TriTriIntersectTest(const Vec3d& v0, const Vec3d& v1, const Vec3d& v2, const Vec3d& nv,
    //    const Vec3d& u0, const Vec3d& u1, const Vec3d& u2, const Vec3d& nu,
    //    int& startType, int& endType, Vec3d& start, Vec3d& end);


}