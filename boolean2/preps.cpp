#include "precompile.h"
#include "preps.h"
#include "xmemory.h"
#include "adaptive.h"

namespace Boolean
{

//#define USE_CGAL_PREDICATES
#define USE_CGAL_PREDICATES_CHECK

#ifdef USE_CGAL_PREDICATES_CHECK
#undef USE_CGAL_PREDICATES
#endif

#ifdef XR_PROFILE

#ifdef XR_DEBUG
#undef XR_DEBUG
#endif

#ifdef USE_CGAL_PREDICATES
#undef USE_CGAL_PREDICATES
#endif

#ifdef USE_CGAL_PREDICATES_CHECK
#undef USE_CGAL_PREDICATES_CHECK
#endif

#endif

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

    /// ***************************************************/
    /// XPlane

    void XPlane::construct_from_three_vertices(const cyPointT &p, const cyPointT &q, const cyPointT & r)
    {
        construct_from_one_vertex_two_edges(p, q - p, r - p);

        assert(has_on(p));
        assert(has_on(q));
        assert(has_on(r));
    }

    void XPlane::construct_from_one_vertex_two_edges(const cyPointT & p, const cyPointT & e0, const cyPointT & e1)
    {
        assert(fp_filter_check(reinterpret_cast<const Real*>(&e0), FP_EDGE_CHECK));
        assert(fp_filter_check(reinterpret_cast<const Real*>(&e1), FP_EDGE_CHECK));
        assert(fp_filter_check(reinterpret_cast<const Real*>(&p), FP_FACTOR));

        XPlaneBase *base = register_base();

        cyPointT* thiz = reinterpret_cast<cyPointT*>(base->data_);
        *thiz = (e0).Cross(e1);
        base->data_[3] = -thiz->Dot(p);

        assert(has_on(p));
    }

    void XPlane::construct_coicident_plane(const PlaneLine & l, const cyPointT & p)
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

        XPlaneBase* base = register_base();
        for (int i = 0; i < 3; i++)
        {
            base->at(i) = 0;
        }

        base->at(maxIndex) = std::copysign(1.0, approxNormal[maxIndex]);
        base->at(3) = base->at(maxIndex) > 0 ? -p[maxIndex] : p[maxIndex];
        assert(has_on(p));
    }

    bool XPlane::parallel(const XPlane & p) const
    {
        const Real* data[2] = { get_data(), p.get_data() };
        //const Real* mat[2] = { data_, p.data_ };
        if (mat2x2det(data) != 0) return false;

        const Real* mat2[2] = { data[0] + 1, data[1] + 1 };
        if (mat2x2det(mat2) != 0) return false;

        const Real mat3k[2][2] = {
            { data[0][2], data[0][0] },
            { data[1][2], data[1][0] }
        };
        const Real* mat3[2] = { mat3k[0], mat3k[1] };
        if (mat2x2det(mat3) != 0) return false;

        return true;
    }

    bool XPlane::normal_equals(const XPlane & p) const
    {
        const Real* data[2] = { get_data(), p.get_data() };
        if (data[0][0] * data[1][0] < 0 ||
            data[0][1] * data[1][1] < 0 ||
            data[0][2] * data[1][2] < 0)
            return false;

        return parallel(p);
    }

    const XPlaneBase & XPlane::get_base() const
    {
        assert(is_valid());
        return xplane(std::abs(id_) - 1);
    }

    XPlaneBase* XPlane::register_base()
    {
        XPlaneBase *plane_base;
        id_ = assign_new_plane(&plane_base)+1;
#ifdef PREP_DEBUG_INFO
        debug_data_ = reinterpret_cast<decltype(debug_data_)>(plane_base->data());
#endif
        return plane_base;
    }

    Oriented_side XPlane::orientation(const PlanePoint & p) const
    {
        Real res = Boolean::orientation(p.plane(0), p.plane(1), p.plane(2), *this);
        //const Real* mat[4] = { p.plane(0).get_data(), p.plane(1).get_data(),
        //    p.plane(2).get_data(), get_data() };
        //Real res = mat4x4det(mat);
        //res *= p.plane(0).signd() * p.plane(1).signd() *
        //    p.plane(2).signd() * signd();
        if (res > 0) return ON_POSITIVE_SIDE;
        if (res < 0) return ON_NEGATIVE_SIDE;
        return ON_ORIENTED_BOUNDARY;
    }

    Oriented_side XPlane::orientation(const cyPointT &p) const
    {
        Oriented_side side = OS_WRONG;

        assert(fp_filter_check(reinterpret_cast<const Real*>(&p), FP_FACTOR));

        const Real* data = get_data();
        const cyPointT* thiz = reinterpret_cast<const cyPointT*>(data);
        double res = thiz->Dot(p);
        if (res > -data[3])
        {
            side = ON_POSITIVE_SIDE;
        }
        else if (res < -data[3])
        {
            side = ON_NEGATIVE_SIDE;
        }
        else
        {
            side = ON_ORIENTED_BOUNDARY;
        }

        if (is_inverse())
        {
            switch (side)
            {
            case ON_NEGATIVE_SIDE:
                return ON_POSITIVE_SIDE;
            case ON_ORIENTED_BOUNDARY:
                return ON_ORIENTED_BOUNDARY;
            case ON_POSITIVE_SIDE:
                return ON_NEGATIVE_SIDE;
            default:
                throw 1;
            }
        }
        else
        {
            return side;
        }
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

    /// ***************************************************/
    /// PlaneLine

    int PlaneLine::linear_order_unsafe(const XPlane & a, const XPlane & b) const
    {
        //const Real* mat[4] = { planes_[0].get_data(),
        //    planes_[1].get_data(), b.get_data(), a.get_data() };

        //Real res = mat4x4det(mat);
        //res *= planes_[0].signd() * planes_[1].signd() * a.signd() * b.signd();
        assert(dot(a) > 0 && dot(b) > 0);

        Real res = Boolean::orientation(planes_[0], planes_[1], b, a);
        if (res > 0) return 1;
        if (res < 0) return -1;
        return 0;
    }

    int PlaneLine::linear_order_unsafe(const XPlane & a, const PlanePoint & b) const
    {
        assert(dot(a) > 0 && has_on(b));

        Oriented_side res = a.orientation(b);
        switch (res)
        {
        case ON_NEGATIVE_SIDE:
            return -1;
        case ON_ORIENTED_BOUNDARY:
            return 0;
        case ON_POSITIVE_SIDE:
            return 1;
        default:
            throw 1;
        }
    }

    int PlaneLine::linear_order_unsafe(const PlanePoint & a, const XPlane & b) const
    {
        assert(has_on(a) && dot(b) > 0);

        Oriented_side res = b.orientation(a);
        switch (res)
        {
        case ON_NEGATIVE_SIDE:
            return 1;
        case ON_ORIENTED_BOUNDARY:
            return 0;
        case ON_POSITIVE_SIDE:
            return -1;
        default:
            throw 1;
        }
    }


    int PlaneLine::linear_order_unsafe(const XPlane & a, const cyPointT & b) const
    {
        assert(dot(a) > 0 && has_on(b));

        Oriented_side res = a.orientation(b);
        switch (res)
        {
        case ON_NEGATIVE_SIDE:
            return -1;
        case ON_ORIENTED_BOUNDARY:
            return 0;
        case ON_POSITIVE_SIDE:
            return 1;
        default:
            throw 1;
        }
    }

    int PlaneLine::linear_order_unsafe(const cyPointT & a, const XPlane & b) const
    {
        assert(has_on(a) && dot(b) > 0);

        Oriented_side res = b.orientation(a);
        switch (res)
        {
        case ON_NEGATIVE_SIDE:
            return 1;
        case ON_ORIENTED_BOUNDARY:
            return 0;
        case ON_POSITIVE_SIDE:
            return -1;
        default:
            throw 1;
        }
    }

    bool PlaneLine::linear_coincident_no_check(const XPlane & a, const XPlane & b) const
    {
        assert(dot(a) != 0 && dot(b) != 0);
        return orientation(planes_[0], planes_[1], a, b) == 0;
    }

    int PlaneLine::linear_order(const PlanePoint & a, const PlanePoint & b) const
    {
        //for (int i = 0; i < 3; i++)
        //{
        //    if (orientation(planes_[0], planes_[1], a.plane(i)) != 0)
        //    {
        //        pa = a.plane(i);
        //        break;
        //    }
        //}

        //for (int i = 0; i < 3; i++)
        //{
        //    if (orientation(planes_[0], planes_[1], b.plane(i)) != 0)
        //    {
        //        pb = b.plane(i);
        //        break;
        //    }
        //}

        XPlane pa, pb;
        pa = pick_positive_vertical_plane(a);
        pb = pick_positive_vertical_plane(b);
        assert(pa.is_valid() && pb.is_valid());

        return linear_order_unsafe(pa, pb);
    }

    int PlaneLine::linear_order(const PlanePoint & a, const XPlane & b) const
    {
        XPlane b_copy = b;
        make_positive(b_copy);
        return linear_order_unsafe(a, b_copy);
    }

    int PlaneLine::linear_order(const XPlane & a, const PlanePoint & b) const
    {
        XPlane a_copy = a;
        make_positive(a_copy);
        return linear_order_unsafe(a_copy, b);
    }

    int PlaneLine::linear_order(const cyPointT& a, const cyPointT& b) const
    {
        assert(has_on(a) & has_on(b));

        cyPointT exact_vec = b - a;
        if (exact_vec.LengthSquared() == Real(0)) return 0;


        const Real* mat[3] = { planes_[0].get_data(), planes_[1].get_data(), (Real*)&exact_vec };

        Real res = mat3x3det(mat);
        res *= planes_[0].signd() * planes_[1].signd();
        if (res > 0) return 1;
        else return -1;
    }

    int PlaneLine::linear_order(const PlanePoint & a, const cyPointT & b) const
    {
        XPlane pa = pick_positive_vertical_plane(a);
        return linear_order_unsafe(pa, b);
    }

    int PlaneLine::linear_order(const cyPointT & a, const PlanePoint & b) const
    {
        XPlane pb = pick_positive_vertical_plane(b);
        return linear_order_unsafe(a, pb);
    }

    int PlaneLine::linear_order(const XPlane & a, const cyPointT & b) const
    {
        XPlane a2 = a;
        make_positive(a2);

        return linear_order_unsafe(a2, b);
    }

    int PlaneLine::linear_order(const cyPointT & a, const XPlane & b) const
    {
        XPlane b2 = b;
        make_positive(b2);

        return linear_order_unsafe(a, b2);
    }

    int PlaneLine::linear_order(const XPlane & a, const XPlane & b) const 
    {
        XPlane a2 = a, b2 = b;
        make_positive(a2);
        make_positive(b2);

        return linear_order_unsafe(a2, b2);
    }

    void PlaneLine::make_positive(XPlane & input) const
	{
        Real res = dot(input);
        if (res < 0) input.inverse();
	}

    Real PlaneLine::dot(const XPlane &input) const
    {
        return orientation(planes_[0], planes_[1], input);
    }

    Real PlaneLine::dot(const Real *p) const
    {
        const Real* mat[3] = { 
            planes_[0].get_data(), 
            planes_[1].get_data(), 
            p };
        Real res = mat3x3det(mat);

        res *= planes_[0].signd() * planes_[1].signd();
        return res;
    }

    XPlane PlaneLine::pick_positive_vertical_plane(const cyPointT & p) const
    {
        XPlane result;
        result.construct_coicident_plane(*this, p);
        return result;
    }

    XPlane PlaneLine::pick_positive_vertical_plane(const PlanePoint & p) const
    {
        for (int j = 0; j < 3; j++)
        {
            Real fres = dot(p.plane(j));
            if (fres == Real(0)) continue;

            XPlane result;
            if (fres > 0)
            {
                result = p.plane(j);
                assert(dot(result) > 0);
                return result;
            }
            else if (fres < 0)
            {
                result = p.plane(j).opposite();
                assert(dot(result) > 0);
                return result;
            }
        }
        throw std::exception("cannnot find a proper plane");
    }

    void PlaneLine::inverse()
    {
        planes_[1].inverse();
    }

    cyPointT PlaneLine::approxNormal() const
    {
        vec3 normal;
        vec3_mul_cross(normal, planes_[0].get_data(), planes_[1].get_data());
        vec3_scale(normal, normal, (planes_[0].signd() * planes_[1].signd()));
        return cyPointT(normal);
    }

    /// ***************************************************/
    /// PlanePoint

    PlanePoint::PlanePoint(const XPlane & a, const XPlane & b, const XPlane & c) :
        planes_{ a, b, c }
    {
        assert(check_positive());
#ifdef PREP_DEBUG_INFO
        auto res = convertToPoint<Depick>(a, b, c);
        coord[0] = res.x();
        coord[1] = res.y();
        coord[2] = res.z();
#endif
    }

    bool PlanePoint::value_equals(const PlanePoint &p) const
    {
        return plane(0).has_on(p)
            && plane(1).has_on(p)
            && plane(2).has_on(p);
    }

    bool PlanePoint::value_equals(const cyPointT & p) const
    {
        return plane(0).has_on(p)
            && plane(1).has_on(p)
            && plane(2).has_on(p);
    }

    bool PlanePoint::check_positive() const
    {
        return orientation(planes_[0], planes_[1], planes_[2]) > 0;
    }

    /// ***************************************************/
    /// other

    Real orientation(const XPlane & p, const XPlane & q, const XPlane & input)
    {
        const Real* mat[3] = { p.get_data(), q.get_data(), input.get_data() };
        Real res = mat3x3det(mat);

        res *= p.signd() * q.signd() * input.signd();
        return res;
    }

    Real orientation(const XPlane & p, const XPlane & q, const XPlane & r, const XPlane & s)
    {
        const Real* mat[4] = { p.get_data(), q.get_data(), r.get_data(), s.get_data() };
        Real res = mat4x4det(mat);

        res *= p.signd() * q.signd() * r.signd() * s.signd();
        return res;
    }
}