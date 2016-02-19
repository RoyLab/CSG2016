#pragma once
#include <CGAL\Bbox_3.h>
#include <CGAL\Iso_cuboid_3.h>
#include <CGAL\Triangle_3.h>


namespace myext
{
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

    inline CGAL::Bbox_3 enlarge(const CGAL::Bbox_3& bbox, double factor)
    {
        return CGAL::Bbox_3(bbox.xmin() - factor, bbox.ymin() - factor,
            bbox.zmin() - factor, bbox.xmax() - factor,
            bbox.ymax() - factor, bbox.zmax() - factor);
    }

    template <class _R>
    inline bool is_inside_box(const CGAL::Iso_cuboid_3<_R>& boxBig, const CGAL::Iso_cuboid_3<_R>& boxSmall)
    {
        return (boxBig.has_on_bounded_side(boxSmall.min()) &&
            boxBig.has_on_bounded_side(boxSmall.max()));
    }

    template <class _R>
    inline typename _R::FT dot(const CGAL::Vector_3<_R>& v, const CGAL::Point_3<_R>& p)
    {
        return v.x() * p.x() + v.y() * p.y() + v.z() * p.z();
    }

    template <class _R>
    inline typename _R::FT dot(const CGAL::Point_3<_R>& p, const CGAL::Vector_3<_R>& v)
    {
        return v.x() * p.x() + v.y() * p.y() + v.z() * p.z();
    }

    template <class _R>
    inline int max_coord(const CGAL::Vector_3<_R>& v)
    {
        typedef typename _R::FT FT;
        int index = 0;

        FT max = fabs(LineDir[0]);
        FT b = fabs(LineDir[1]);
        FT c = fabs(LineDir[2]);

        if (b > max) { max = b; index = 1; }
        if (c > max) { max = c; index = 2; }

        return index;
    }

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
        typedef typename _R::FT FT;

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

    //int TriTriIntersectTest(const Vec3d& v0, const Vec3d& v1, const Vec3d& v2, const Vec3d& nv,
    //    const Vec3d& u0, const Vec3d& u1, const Vec3d& u2, const Vec3d& nu,
    //    int& startType, int& endType, Vec3d& start, Vec3d& end);

}