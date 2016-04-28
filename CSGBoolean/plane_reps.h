#pragma once
#include "CGALext.h"
#include "macroutil.h"
#include "adaptive.h"
#include "csgdefs.h"
#include <cstring>
#include <cmath>
#include <CGAL\intersections.h>

namespace CSG
{
    const double P216 = std::pow(2.0, 16);
    const double P217 = std::pow(2.0, 17);
    const double P218 = std::pow(2.0, 18);

    extern double FP_FACTOR;


    template <class Vector>
    inline Vector filter(Vector& p, double factor)
    {
        return Vector(GS::fp_filter(p.x(), factor),
            GS::fp_filter(p.y(), factor), GS::fp_filter(p.z(), factor));
    }

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
        PR_None,
        PR_In,
        PR_Out,
        PR_OnBoundary
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
            CGAL::Plane_3<_R>(p, CGAL::cross_product(q - p, r - p)){}

        Plane_ext(const Point& p, const Vector& e0, const Vector& e1):
            CGAL::Plane_3<_R>(p, CGAL::cross_product(e0, e1)){}

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
            ReportError("");
            return 0;
        }
    };

    template <class _R>
    struct PBPoint
    {
        DECLARE_CGAL_KERNEL_CLASS
        typedef Plane_ext<_R> PlaneExt;
        template <class T>
        friend std::ostream& operator<<(std::ostream& out, const PBPoint<T>& p);
        PBPoint(){}
        PBPoint(const Plane& p, const Plane& q, const Plane& r, bool keepPositive = true)
        {
            planes[0] = p;
            planes[1] = q;
            planes[2] = r;

            if (keepPositive)
            {
                if (CGAL::determinant(p.orthogonal_vector(),
                    q.orthogonal_vector(), r.orthogonal_vector()) < 0.0)
                    planes[2] = r.opposite();
            }

            assert(CGAL::determinant(p.orthogonal_vector(),
                q.orthogonal_vector(), r.orthogonal_vector()) > 0.0);

            computeCoord();
        }

        PBPoint(const PBPoint& p)
        {
            for (size_t i = 0; i < 3; i++)
                planes[i] = p.planes[i];
            coord = p.coord;
        }

        bool operator==(const PBPoint& p) const
        {
            double4x4 mat;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 4; j++)
                    mat[i][j] = planes[i][j];

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 4; j++)
                    mat[3][j] = p.planes[i][j];

                if (GS::adaptiveDet4x4Sign(mat) != 0)
                    return false;
            }

            //std::cout << "\nequal event:\n";
            //std::cout << *this;
            //std::cout << p;
            return true;
        }

        double estimateSquaredDistance(const PBPoint& p) const
        {
            return (coord - p.coord).squared_distance();
        }

        RelationToPlane classifyByPlane(const Plane_ext<_R>& p) const
        {
            return p.classifyPointToPlane(planes[0], planes[1], planes[2]);
        }

        void computeCoord()
        {
            auto result = CGAL::intersection(planes[0], planes[1], planes[2]);
            const Point* p = boost::get<Point>(&*result);
            coord = *p;
        }

        const PlaneExt* getPlanes() const { return planes; }
        const PlaneExt& getPlane(size_t i) const { return planes[i]; }
        const Point& getCoord() const { return coord; }

    private:
        PlaneExt planes[3];
        Point coord;
    };

    template <class _R>
    std::ostream& operator<<(std::ostream& out, const PBPoint<_R>& p)
    {
        std::cout << p.planes[0] << std::endl;
        std::cout << p.planes[1] << std::endl;
        std::cout << p.planes[2] << std::endl;
        std::cout << p.coord << std::endl;
        return out;
    }

    /* bounding planes point to the internal of the triangle */
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

            auto pseudo = filter(normal * 0.1, FP_FACTOR);

            m_bps[0] = PlaneExt(q, pseudo, r - q);
            m_bps[1] = PlaneExt(r, pseudo, p - r);
            m_bps[2] = PlaneExt(p, pseudo, q - p);

            assert(m_sp.has_on(p));
            assert(m_sp.has_on(q));
            assert(m_sp.has_on(r));
            assert(m_bps[0].has_on(q));
            assert(m_bps[0].has_on(r));
            assert(m_bps[1].has_on(r));
            assert(m_bps[1].has_on(p));
            assert(m_bps[2].has_on(p));
            assert(m_bps[2].has_on(q));
        }

        PBPoint<_R> point(int idx)
        {
            return PBPoint<_R>(m_bps[(idx + 1) % 3], m_bps[(idx + 2) % 3], m_sp);
        }

        bool insideBPs(const PBPoint<_R>& p) const
        {
            for (int i = 0; i < 3; i++)
            {
                if (p.classifyByPlane(m_bps[i]) == Behind)
                    return false;
            }
            return true;
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
        double4x4 mat{ p.a(), p.b(), p.c(), p.d(),
            q.a(), q.b(), q.c(), q.d(),
            r.a(), r.b(), r.c(), r.d(),
            s.a(), s.b(), s.c(), s.d() };
        double sign1 = GS::adaptiveDet4x4Sign(mat);
#ifdef _DEBUG
        CGAL::Sign sign2 = CGAL::sign_of_determinant(p.a(), p.b(), p.c(), p.d(),
            q.a(), q.b(), q.c(), q.d(),
            r.a(), r.b(), r.c(), r.d(),
            s.a(), s.b(), s.c(), s.d());

        double d = std::abs(CGAL::determinant(p.a(), p.b(), p.c(), p.d(),
            q.a(), q.b(), q.c(), q.d(),
            r.a(), r.b(), r.c(), r.d(),
            s.a(), s.b(), s.c(), s.d()));
        assert(sign1 == 0.0 && (sign2 == CGAL::ZERO || d < 1e-10) ||
            sign1 < 0.0 && sign2 == CGAL::NEGATIVE ||
            sign1 > 0.0 && sign2 == CGAL::POSITIVE);
#endif
        return sign1;
    }

    template <class Vector>
    double orientation(const Vector& p, const Vector& q, const Vector& r)
    {
        double3x3 mat{ p.x(), p.y(), p.z(),
            q.x(), q.y(), q.z(),
            r.x(), r.y(), r.z() };

        double sign1 = GS::adaptiveDet3x3Sign(mat);
#ifdef _DEBUG
        CGAL::Sign sign2 = CGAL::sign_of_determinant(p.x(), p.y(), p.z(),
            q.x(), q.y(), q.z(),
            r.x(), r.y(), r.z());

        double d = std::abs(CGAL::determinant(p.x(), p.y(), p.z(),
            q.x(), q.y(), q.z(),
            r.x(), r.y(), r.z()));
        assert(sign1 == 0.0 && (sign2 == CGAL::ZERO || d < 1e-10) ||
            sign1 < 0.0 && sign2 == CGAL::NEGATIVE ||
            sign1 > 0.0 && sign2 == CGAL::POSITIVE);
#endif
        return sign1;
    }

    template <class _R>
    static inline bool same_orientation(const PBPoint<_R>& p0, const PBPoint<_R>& p1, const PBPoint<_R>& p2, const CGAL::Vector_3<_R>& normal)
    {
        auto e1 = p1.getCoord() - p0.getCoord();
        auto e2 = p2.getCoord() - p1.getCoord();

        return CGAL::orientation(e1, e2, normal) == CGAL::POSITIVE;
    }

    template <class RT>
    static inline double determinant3x3(const RT& a00, const RT& a01, const RT& a02,
        const RT& a10, const RT& a11, const RT& a12,
        const RT& a20, const RT& a21, const RT& a22)
    {
        double3x3 mat = { a00, a01, a02, a10, a11, a12, a20, a21, a22 };
        double ss = GS::adaptiveDet3x3Sign(mat);
#ifdef _DEBUG
        //CGAL::Sign s = CGAL::sign_of_determinant(a00, a01, a02, a10, a11, a12, a20, a21, a22);
        //assert(ss == 0.0 && s == CGAL::ZERO ||
        //    ss < 0.0 && s == CGAL::NEGATIVE ||
        //    ss > 0.0 && s == CGAL::POSITIVE);
#endif
        return ss;
    }

    template <class RT>
    static inline double determinant4x4(const RT& a00, const RT& a01, const RT& a02, const RT& a03,
        const RT& a10, const RT& a11, const RT& a12, const RT& a13,
        const RT& a20, const RT& a21, const RT& a22, const RT& a23,
        const RT& a30, const RT& a31, const RT& a32, const RT& a33)
    {
        double3x3 mat = { a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33 };
        double ss = GS::adaptiveDet4x4Sign(mat);
#ifdef _DEBUG
        //CGAL::Sign s = CGAL::sign_of_determinant(a00, a01, a02, a03, a10, a11, a12, a13, a20, a21, a22, a23, a30, a31, a32, a33);
        //assert(ss == 0.0 && s == CGAL::ZERO ||
        //    ss < 0.0 && s == CGAL::NEGATIVE ||
        //    ss > 0.0 && s == CGAL::POSITIVE);
#endif
        return ss;
    }

    static inline int determinePhase(double abcos, double absin)
    {
        if (abcos >= 0.0 && absin >= 0.0)
            return 1;

        if (abcos < 0.0 && absin >= 0.0)
            return 2;

        if (abcos < 0.0 && absin < 0.0)
            return 3;

        else return 4;
    }

}