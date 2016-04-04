#pragma once
#include "CGALext.h"
#include "macroutil.h"
#include "adaptive.h"
#include "csgdefs.h"
#include <cstring>
#include <CGAL\intersections.h>

namespace CSG
{
    enum RelationToPlane{
        Front,
        On,
        Behind,
        Straddling
    };

    enum OutputSymbol{
        Empty,
        B,
        HB,
        PB,
        PBB,
    };

    enum PosRelation{
        In,
        Out,
        OnBoundary
    };

    template <class _R>
    class Plane_ext :
        public CGAL::Plane_3<_R>
    {
        DECLARE_CGAL_KERNEL_CLASS
    public:
        Plane_ext(){}
        Plane_ext(const CGAL::Plane_3<_R>& p)
        { memcpy(this, &p, sizeof(Plane_ext)); }

        Plane_ext(const Point& p, const Point& q, const Point& r):
            CGAL::Plane_3<_R>(p, CGAL::cross_product(q - p, r - p))
        {}

        RelationToPlane classifyPointToPlane(const Plane_ext & p, const Plane_ext & q, const Plane_ext & r) const
        {
            double3x3 ptMat;
            for (int i = 0; i < 3; i++)
            {
                ptMat[0][i] = p[i];
                ptMat[1][i] = q[i];
                ptMat[2][i] = r[i];
            }

            double4x4 matPlane;
            for (int i = 0; i < 4; i++)
            {
                matPlane[0][i] = p[i];
                matPlane[1][i] = q[i];
                matPlane[2][i] = r[i];
                matPlane[3][i] = this->operator[](i);
            }

            auto dist = GS::adaptiveDet3x3Sign(ptMat)*GS::adaptiveDet4x4Sign(matPlane);
            if (dist > 0.0) return Front;
            if (dist < 0.0) return Behind;
            return On;
        }

        typename _R::FT operator[](int i) const
        {
            assert(i < 4 && i >= 0);
            switch (i)
            {
            case 0: return a();
            case 1: return b();
            case 2: return c();
            case 3: return d();
            }
            ReportError();
            return 0;
        }
    };

    template <class _R>
    struct PBPoint
    {
        DECLARE_CGAL_KERNEL_CLASS
        typedef Plane_ext<_R> PlaneExt;

        PBPoint(){}
        PBPoint(const Plane& p, const Plane& q, const Plane& r, bool keepPositive = true)
        {
            planes[0] = p;
            planes[1] = q;
            planes[2] = r;

            if (keepPositive)
            {
                if (CGAL::determinant(p.orthogonal_vector(),
                    q.orthogonal_vector(), r.orthogonal_vector()) < 0.f)
                    planes[2] = r.opposite();
            }

            assert(CGAL::determinant(p.orthogonal_vector(),
                q.orthogonal_vector(), r.orthogonal_vector()) > 0.f);

            computeCoord();
        }

        //Point estimate_coords() const
        //{
        //    return CGAL::intersection(planes[0], planes[1], planes[2]);
        //}

        bool operator==(const PBPoint& p) const
        {
            double4x4 mat;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 4; j++)
                    mat[i][j] = p.planes[i][j];

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 4; j++)
                    mat[3][j] = p.planes[i][j];

                if (GS::adaptiveDet4x4Sign(mat) != 0)
                    return false;
            }
            return true;
        }


        void computeCoord()
        {
            auto result = CGAL::intersection(planes[0], planes[1], planes[2]);
            const Point* p = boost::get<Point>(&*result);
            coord == *p;
        }

        const PlaneExt* getPlanes() const { return planes; }
        const Point& getCoord() const { return coord; }

    private:
        PlaneExt planes[3];
        Point coord;
    };

    template <class _R>
    struct PBTriangle
    {
        DECLARE_CGAL_KERNEL_CLASS
        typedef Plane_ext<_R> PlaneExt;

        COMMON_PROPERTY(PlaneExt, sp);
        PlaneExt m_bps[3];

        PBTriangle(){}

        PBTriangle(const Triangle& t)
        {
            new(this)PBTriangle(t.vertex(0), t.vertex(1), t.vertex(2));
        }

        PBTriangle(const Point& p, const Point& q, const Point& r)
        {
            m_sp = PlaneExt(p, q, r);

            Vector normal = m_sp.orthogonal_vector();
            normal = normal / sqrt(normal.squared_length());

            m_bps[0] = PlaneExt(q, r, q + normal);
            m_bps[1] = PlaneExt(r, p, r + normal);
            m_bps[2] = PlaneExt(p, q, p + normal);
        }

        PBPoint<_R> point(int idx)
        {
            return PBPoint<_R>(m_bps[(idx + 1) % 3], m_bps[(idx + 2) % 3], m_sp);
        }
    };


    template <class _R>
    struct PBPolygon{
        //typedef K _R;

        DECLARE_CGAL_KERNEL_CLASS
        typedef Plane_ext<_R> PlaneExt;

        COMMON_PROPERTY(PlaneExt, sp);
        std::vector<PlaneExt> m_bps;

        PBPolygon(){}
        PBPolygon(const PBTriangle<_R>& triangle)
        {
            m_sp = triangle.get_sp();
            for (int i = 0; i < 3; i++)
                m_bps.push_back(triangle.m_bps[i]);
        }

        int Size() const { return m_bps.size(); }
        void Clear() { m_bps.clear(); }
        RelationToPlane ClassifyPloygonToPlane(const PlaneExt& plane) const
        {
            int numInFront = 0;
            int numInBack = 0;

            for (int i = 0; i < m_bps.size(); i++)
            {
                int pBpId = i > 0 ? i - 1 : m_bps.size() - 1;
                RelationToPlane PtToPlane = plane.classifyPointToPlane(m_sp, m_bps[pBpId], m_bps[i]);
                if (PtToPlane == Front)
                    numInFront++;
                else if (PtToPlane == Behind)
                    numInBack++;
            }
            if (numInFront != 0 && numInBack != 0)
                return Straddling;
            if (numInFront != 0)
                return Front;
            if (numInBack != 0)
                return Behind;
            return On;
        }

        RelationToPlane ClassifyEdgeToPlane(const PlaneExt& plane, int idx) const;
        RelationToPlane ClipByPlane(const PlaneExt& bp, PBPolygon& front) const; // only output the polygon in the front of  clipping plane
        RelationToPlane ClipByPlaneNoFront(PlaneExt& bp); // only output the polygon in the front of  clipping plane

        RelationToPlane ClipByPlane(const PlaneExt& bp, PBPolygon& front, PBPolygon& back) const
        {
            if (m_sp == bp)
                return On;

            auto n = m_bps.size();
            RelationToPlane *relation = new RelationToPlane[n];
            for (int i = 0; i < n; i++)
                relation[i] = bp.classifyPointToPlane(m_sp, m_bps[i], m_bps[(i + 1) % n]);

            bool bFrontInserted = false;
            bool bBackInserted = false;

            for (int i = 0; i < n; i++)
            {
                // classify  point Vx-1 , Vx , Vx+1 to Plane  bp
                RelationToPlane  prevPtPos = relation[(i + n - 2) % n];
                RelationToPlane currentPtPos = relation[(i + n - 1) % n];
                RelationToPlane nextPtPos = relation[i];

                OutputSymbol signal = LookupEncodingTable(prevPtPos, currentPtPos, nextPtPos);
                switch (signal)
                {
                case B:
                    front.m_bps.push_back(m_bps[i]);
                    break;
                case HB:
                    if (!bFrontInserted)
                    {
                        front.m_bps.push_back(bp);
                        bFrontInserted = true;
                    }
                    front.m_bps.push_back(m_bps[i]);
                    break;
                default:
                    break;
                }

                // output back plane
                OutputSymbol backSignal = LookupEncodingTable(ReverseRelation(prevPtPos), ReverseRelation(currentPtPos), ReverseRelation(nextPtPos));
                switch (backSignal)
                {
                case B:
                    back.m_bps.push_back(m_bps[i]);
                    break;
                case HB:
                    if (!bBackInserted)
                    {
                        PlaneExt hp(bp.opposite());
                        back.m_bps.push_back(hp);
                        bBackInserted = true;
                    }
                    back.m_bps.push_back(m_bps[i]);
                    break;
                default:
                    break;
                }

            }
            delete[] relation;
            auto res = On;
            if (front.m_bps.size())
            {
                front.m_sp = m_sp;
                res = Front;
            }
            if (back.m_bps.size())
            {
                back.m_sp = m_sp;
                res = Behind;
            }
            if (front.m_bps.size() && back.m_bps.size())
                return Straddling;

            return res;
        }

        RelationToPlane ClipByPlaneNoFront(Plane& bp, PBPolygon& back);
        void Negate();

    private:
        OutputSymbol  LookupEncodingTable(RelationToPlane prevPtPos, RelationToPlane currentPtPos, RelationToPlane nextPtPos) const
        {
            if (currentPtPos == Front)
                return B;
            if (currentPtPos == Behind)
            {
                if (nextPtPos == Front)
                    return HB;
                return Empty;
            }
            if (nextPtPos == Front)
            {
                if (prevPtPos == Front)
                    return B;
                return  HB;
            }
            else
                return Empty;
        }

        RelationToPlane ReverseRelation(RelationToPlane relation) const
        {
            switch (relation)
            {
            case Front:
                return Behind;
            case Behind:
                return Front;
            default:
                return relation;
            }
        }
    };


    template <class Plane>
    double orientation(const Plane& p, const Plane& q, const Plane& r, const Plane& s)
    {
        return CGAL::determinant(p.a(), p.b(), p.c(), p.d(),
            q.a(), q.b(), q.c(), q.d(),
            r.a(), r.b(), r.c(), r.d(),
            s.a(), s.b(), s.c(), s.d());
    }
}