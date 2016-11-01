#include "precompile.h"
#include "preps.h"
#include "xmemory.h"
#include "adaptive.h"

#define USE_CGAL_PREDICATES

namespace Boolean
{
    XPlane::XPlane(const cyPointT & p, const cyPointT & q, const cyPointT & r)
    {
        xplanes().emplace_back(p, q, r);
        setId(xplanes().size()-1);
        assert(has_on(p));
        assert(has_on(q));
        assert(has_on(r));
#ifdef PREP_DEBUG_INFO
        debug();
#endif
    }

    const XPlaneBase & XPlane::base() const
	{
		return xplane(std::abs(id) - 1);
	}

    Oriented_side XPlane::orientation(const XPoint & p) const
    {
        const Real* mat[4] = { p.plane(0).data(), p.plane(1).data(),
            p.plane(2).data(), data() };
#ifdef USE_CGAL_PREDICATES
        Real res = cgalExact4x4(mat);
#else
        Real res = adaptiveDet4x4Sign(mat);
#endif
        res *= p.plane(0).signd() * p.plane(1).signd() *
            p.plane(2).signd() * signd();
        if (res > 0) return ON_POSITIVE_SIDE;
        if (res < 0) return ON_NEGATIVE_SIDE;
        return ON_ORIENTED_BOUNDARY;
    }

    Oriented_side XPlane::orientation(const cyPointT &p) const
    {
        assert(fp_filter_check(reinterpret_cast<const Real*>(&p), FP_FACTOR));

        const cyPointT* thiz = reinterpret_cast<const cyPointT*>(data());
        double res = thiz->Dot(p);
        res *= signd();
        if (res > -data()[3] * signd()) return ON_POSITIVE_SIDE;
        if (res < -data()[3] * signd()) return ON_NEGATIVE_SIDE;
        return ON_ORIENTED_BOUNDARY;
    }

    void XPlane::setFromPEE(const cyPointT & p, const cyPointT & e0, const cyPointT & e1)
    {
        xplanes().emplace_back();
        xplanes().back().setFromPEE(p, e0, e1);
        setId(xplanes().size() - 1);
#ifdef PREP_DEBUG_INFO
        debug();
#endif
    }

    int XLine::linearOrderNoCheck(const XPlane & a, const XPlane & b) const
    {
        const Real* mat[4] = { m_planes[0].data(),
            m_planes[1].data(), b.data(), a.data() };

#ifdef USE_CGAL_PREDICATES
        Real res = cgalExact4x4(mat);
#else
        Real res = adaptiveDet4x4Sign(mat);
#endif
        res *= m_planes[0].signd() * m_planes[1].signd() * a.signd() * b.signd();
        if (res > 0) return 1;
        if (res < 0) return -1;
        return 0;
    }

    int XLine::linearOrder(const XPoint & a, const XPoint & b) const
    {
        XPlane pa, pb;
        for (int i = 0; i < 3; i++)
        {
            if (sign(m_planes[0], m_planes[1], a.plane(i)) != 0)
            {
                pa = a.plane(i);
                break;
            }
        }

        for (int i = 0; i < 3; i++)
        {
            if (sign(m_planes[0], m_planes[1], b.plane(i)) != 0)
            {
                pb = b.plane(i);
                break;
            }
        }
        assert(pa.isValid() && pb.isValid());
        return linearOrder(pa, pb);
    }

    int XLine::linearOrder(const cyPointT& a, const cyPointT& b) const
    {
        cyPointT vec = b - a;
        if (vec.LengthSquared() == Real(0)) return 0;

        const Real* mat[3] = { m_planes[0].data(), m_planes[1].data(), (Real*)&vec };
#ifdef USE_CGAL_PREDICATES
        Real res = cgalExact3x3(mat);
#else
        Real res = adaptiveDet3x3Sign(mat);
#endif
        res *= m_planes[0].signd() * m_planes[1].signd();
        if (res > 0) return 1;
        else return -1;
    }

    int XLine::linearOrder(const XPlane & a, const XPlane & b) const 
    {
        auto a2 = a, b2 = b;
        makePositive(a2);
        makePositive(b2);

        const Real* mat[4] = { m_planes[0].data(),
            m_planes[1].data(), b2.data(), a2.data() };

#ifdef USE_CGAL_PREDICATES
        Real res = cgalExact4x4(mat);
#else
        Real res = adaptiveDet4x4Sign(mat);
#endif
        res *= m_planes[0].signd() * m_planes[1].signd() *
            a2.signd() * b2.signd();

        if (res > 0) return 1;
        if (res < 0) return -1;
        return 0;
    }

    void XLine::makePositive(XPlane & input) const
	{
        Real res = dot(input);
        if (res < 0) input.inverse();
	}

    Real XLine::dot(const XPlane &input) const
    {
        return sign(m_planes[0], m_planes[1], input);
    }

    XPlane XLine::pickPositiveVertical(const XPoint & p) const
    {
        for (int j = 0; j < 3; j++)
        {
            Real fres = dot(p.plane(j));
            if (fres == Real(0)) continue;

            if (fres > 0)
                return p.plane(j);
            else if (fres < 0)
                return p.plane(j).opposite();
        }
        throw std::exception();
        return XPlane();
    }

    XPlaneBase::XPlaneBase(const cyPointT &p, const cyPointT &q, const cyPointT & r)
    {
        setFromPEE(p, q - p, r - p);
    }

    void XPlaneBase::setFromPEE(const cyPointT & p, const cyPointT & e0, const cyPointT & e1)
    {
        assert(fp_filter_check(reinterpret_cast<const Real*>(&e0), FP_EDGE_CHECK));
        assert(fp_filter_check(reinterpret_cast<const Real*>(&e1), FP_EDGE_CHECK));

        cyPointT* thiz = reinterpret_cast<cyPointT*>(this);
        *thiz = (e0).Cross(e1);
        m_data[3] = -thiz->Dot(p);
    }

    bool XPoint::operator==(const XPoint &p) const
    {
        return plane(0).orientation(p) == ON_ORIENTED_BOUNDARY
            && plane(1).orientation(p) == ON_ORIENTED_BOUNDARY
            && plane(2).orientation(p) == ON_ORIENTED_BOUNDARY;
    }

    bool XPoint::operator==(const cyPointT & p) const
    {
        return plane(0).orientation(p) == ON_ORIENTED_BOUNDARY
            && plane(1).orientation(p) == ON_ORIENTED_BOUNDARY
            && plane(2).orientation(p) == ON_ORIENTED_BOUNDARY;

    }
    Real sign(const XPlane & p, const XPlane & q, const XPlane & input)
    {
        const Real* mat[3] = { p.data(), q.data(), input.data() };
#ifdef USE_CGAL_PREDICATES
        Real res = cgalExact3x3(mat);
#else
        Real res = adaptiveDet3x3Sign(mat);
#endif
        res *= p.signd() * q.signd() * input.signd();
        return res;
    }
}