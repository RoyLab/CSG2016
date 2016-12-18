#include "precompile.h"
#include "preps.h"
#include "xmemory.h"
#include "adaptive.h"

namespace Boolean
{
    namespace
    {
        Real mat2x2det(const Real* mat[2])
        {
#ifdef USE_CGAL_PREDICATES
            Real res = cgalExact2x2(mat);
#else
            Real res = adaptiveDet2x2Sign(mat);
#endif
#ifdef USE_CGAL_PREDICATES_CHECK
            Real check = cgalExact2x2(mat);
            assert(check == res || check * res > 0.0);
#endif
            return res;
        }

        Real mat3x3det(const Real* mat[3])
        {
#ifdef USE_CGAL_PREDICATES
            Real res = cgalExact3x3(mat);
#else
            Real res = adaptiveDet3x3Sign(mat);
#endif
#ifdef USE_CGAL_PREDICATES_CHECK
            Real check = cgalExact3x3(mat);
            assert(check == res || check * res > 0.0);
#endif
            return res;
        }

        Real mat4x4det(const Real* mat[4])
        {
#ifdef USE_CGAL_PREDICATES
            Real res = cgalExact4x4(mat);
#else
            Real res = adaptiveDet4x4Sign(mat);
#endif
#ifdef USE_CGAL_PREDICATES_CHECK
            Real check = cgalExact4x4(mat);
            assert(check == res || check * res > 0.0);
#endif
            return res;
        }
    }

//    XPlane::XPlane(const cyPointT & p, const cyPointT & q, const cyPointT & r)
//    {
//        xplanes().emplace_back(p, q, r);
//        set_positive_from_id(xplanes().size()-1);
//#ifdef PREP_DEBUG_INFO
//        debug();
//#endif
//        assert(has_on(p));
//        assert(has_on(q));
//        assert(has_on(r));
//    }
//
//    XPlane::XPlane(const PlaneLine & l, const cyPointT& p)
//    {
//        xplanes().emplace_back(l, p);
//        set_positive_from_id(xplanes().size() - 1);
//#ifdef PREP_DEBUG_INFO
//        debug();
//#endif
//        assert(has_on(p));
//    }

    void XPlane::setBase(const XPlaneBase & base)
    {
        xplanes().push_back(base);
        set_positive_from_id(xplanes().size() - 1);
#ifdef PREP_DEBUG_INFO
        debug();
#endif
    }

    const XPlaneBase & XPlane::base() const
	{
        assert(id_ != 0);
		return xplane(std::abs(id_) - 1);
	}

    Oriented_side XPlane::orientation(const PlanePoint & p) const
    {
        const Real* mat[4] = { p.plane(0).data(), p.plane(1).data(),
            p.plane(2).data(), data() };

        Real res = mat4x4det(mat);

        res *= p.plane(0).signd() * p.plane(1).signd() *
            p.plane(2).signd() * signd();
        if (res > 0) return ON_POSITIVE_SIDE;
        if (res < 0) return ON_NEGATIVE_SIDE;
        return ON_ORIENTED_BOUNDARY;
    }

    Oriented_side XPlane::orientation(const cyPointT &p) const
    {
        Oriented_side res = base().orientation(p);
        if (is_inverse())
        {
            switch (res)
            {
            case ON_NEGATIVE_SIDE:
                return ON_POSITIVE_SIDE;
            case ON_ORIENTED_BOUNDARY:
                return ON_ORIENTED_BOUNDARY;
            case ON_POSITIVE_SIDE:
                return ON_NEGATIVE_SIDE;
            default:
                return OS_WRONG;
            }
        }
        else
            return res;
    }

//    void XPlane::setFromPEE(const cyPointT & p, const cyPointT & e0, const cyPointT & e1)
//    {
//        xplanes().emplace_back();
//        xplanes().back().setFromPEE(p, e0, e1);
//        set_positive_from_id(xplanes().size() - 1);
//#ifdef PREP_DEBUG_INFO
//        debug();
//#endif
//    }

    int PlaneLine::linearOrderNoCheck(const XPlane & a, const XPlane & b) const
    {
        const Real* mat[4] = { planes_[0].data(),
            planes_[1].data(), b.data(), a.data() };

        Real res = mat4x4det(mat);
        res *= planes_[0].signd() * planes_[1].signd() * a.signd() * b.signd();
        if (res > 0) return 1;
        if (res < 0) return -1;
        return 0;
    }

    int PlaneLine::linearOrder(const PlanePoint & a, const PlanePoint & b) const
    {
        XPlane pa, pb;
        for (int i = 0; i < 3; i++)
        {
            if (sign(planes_[0], planes_[1], a.plane(i)) != 0)
            {
                pa = a.plane(i);
                break;
            }
        }

        for (int i = 0; i < 3; i++)
        {
            if (sign(planes_[0], planes_[1], b.plane(i)) != 0)
            {
                pb = b.plane(i);
                break;
            }
        }
        assert(pa.is_valid() && pb.is_valid());
        return linearOrder(pa, pb);
    }

    int PlaneLine::linearOrder(const cyPointT& a, const cyPointT& b) const
    {
        cyPointT vec = b - a;
        if (vec.LengthSquared() == Real(0)) return 0;

        const Real* mat[3] = { planes_[0].data(), planes_[1].data(), (Real*)&vec };

        Real res = mat3x3det(mat);
        res *= planes_[0].signd() * planes_[1].signd();
        if (res > 0) return 1;
        else return -1;
    }

    int PlaneLine::linearOrder(const XPlane & a, const XPlane & b) const 
    {
        auto a2 = a, b2 = b;
        makePositive(a2);
        makePositive(b2);

        const Real* mat[4] = { planes_[0].data(),
            planes_[1].data(), b2.data(), a2.data() };

        Real res = mat4x4det(mat);
        res *= planes_[0].signd() * planes_[1].signd() *
            a2.signd() * b2.signd();

        if (res > 0) return 1;
        if (res < 0) return -1;
        return 0;
    }

    void PlaneLine::makePositive(XPlane & input) const
	{
        Real res = dot(input);
        if (res < 0) input.inverse();
	}

    Real PlaneLine::dot(const XPlane &input) const
    {
        return sign(planes_[0], planes_[1], input);
    }

    XPlane PlaneLine::pickPositiveVertical(const cyPointT & p) const
    {
        return XPlane(XPlaneBase(*this, p));
    }

    XPlane PlaneLine::pickPositiveVertical(const PlanePoint & p) const
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
        throw std::exception("cannnot find a proper plane");
    }

    cyPointT PlaneLine::approxNormal() const
    {
        vec3 normal;
        vec3_mul_cross(normal, planes_[0].data(), planes_[1].data());
        return cyPointT(normal);
    }

    XPlaneBase::XPlaneBase(const PlaneLine & l, const cyPointT & p)
    {
        cyPointT approxNormal = l.approxNormal();
        int maxIndex = 0;
        Real maxVal = std::abs(approxNormal.x);

        for (int i = 1; i < 3; i++)
        {
            if (std::abs(approxNormal.y) > maxVal)
            {
                maxIndex = i;
                maxVal = std::abs(approxNormal[i]);
            }
        }

        for (int i = 0; i < 3; i++) data_[i] = 0;
        data_[maxIndex] = std::copysign(1.0, approxNormal[maxIndex]);
        data_[3] = data_[maxIndex] > 0? -p[maxIndex]: p[maxIndex];
    }

    XPlaneBase::XPlaneBase(const cyPointT &p, const cyPointT &q, const cyPointT & r)
    {
        cyPointT e0 = q - p, e1 = r - p;
        assert(fp_filter_check(reinterpret_cast<const Real*>(&e0), FP_EDGE_CHECK));
        assert(fp_filter_check(reinterpret_cast<const Real*>(&e1), FP_EDGE_CHECK));

        cyPointT* thiz = reinterpret_cast<cyPointT*>(this);
        *thiz = (e0).Cross(e1);
        data_[3] = -thiz->Dot(p);
    }

    //void XPlaneBase::setFromPEE(const cyPointT & p, const cyPointT & e0, const cyPointT & e1)
    //{
    //    assert(fp_filter_check(reinterpret_cast<const Real*>(&e0), FP_EDGE_CHECK));
    //    assert(fp_filter_check(reinterpret_cast<const Real*>(&e1), FP_EDGE_CHECK));

    //    cyPointT* thiz = reinterpret_cast<cyPointT*>(this);
    //    *thiz = (e0).Cross(e1);
    //    m_data[3] = -thiz->Dot(p);
    //}

    Oriented_side XPlaneBase::orientation(const cyPointT & p) const
    {
        assert(fp_filter_check(reinterpret_cast<const Real*>(&p), FP_FACTOR));
        const cyPointT* thiz = reinterpret_cast<const cyPointT*>(data());
        double res = thiz->Dot(p);
        if (res > -data()[3]) return ON_POSITIVE_SIDE;
        if (res < -data()[3]) return ON_NEGATIVE_SIDE;
        return ON_ORIENTED_BOUNDARY;
    }

    bool XPlaneBase::coplanar(const XPlaneBase & p) const
    {
        const Real* mat[2] = { data_, p.data_ };
        if (mat2x2det(mat) != 0) return false;

        const Real* mat2[2] = { data_+1, p.data_+1 };
        if (mat2x2det(mat2) != 0) return false;

        return true;
    }

    bool PlanePoint::operator==(const PlanePoint &p) const
    {
        //return m_pos == p.m_pos;
        return plane(0).orientation(p) == ON_ORIENTED_BOUNDARY
            && plane(1).orientation(p) == ON_ORIENTED_BOUNDARY
            && plane(2).orientation(p) == ON_ORIENTED_BOUNDARY;
    }

    bool PlanePoint::operator==(const cyPointT & p) const
    {
        return plane(0).orientation(p) == ON_ORIENTED_BOUNDARY
            && plane(1).orientation(p) == ON_ORIENTED_BOUNDARY
            && plane(2).orientation(p) == ON_ORIENTED_BOUNDARY;

    }
    Real sign(const XPlane & p, const XPlane & q, const XPlane & input)
    {
        const Real* mat[3] = { p.data(), q.data(), input.data() };
        Real res = mat3x3det(mat);

        res *= p.signd() * q.signd() * input.signd();
        return res;
    }
}