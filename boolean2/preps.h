#pragma once
#ifdef PREP_DEBUG_INFO
#include <CGAL/intersection_3.h>
#endif
#include "global.h"

namespace Boolean
{
	class XPlaneBase
	{
	public:
        XPlaneBase() {}
        XPlaneBase(const PlaneLine& l, const cyPointT& p); // give a p-rep normal and a point, construct a approx plane on p
        XPlaneBase(const cyPointT& p, const cyPointT& q, const cyPointT& r);
        XPlaneBase(const cyPointT& p, const cyPointT& e0, const cyPointT& e1);

        Oriented_side orientation(const cyPointT&) const;
        bool coplanar(const XPlaneBase&) const;

		const Real& a() const { return data_[0]; }
		const Real& b() const { return data_[1]; }
		const Real& c() const { return data_[2]; }
		const Real& d() const { return data_[3]; }

        const Real* data() const { return data_; }
        Real* data() { return data_; }

	protected:
		Real data_[4];
	};

	class XPlane
	{
	public:
		XPlane() : id_(0) {
#ifdef PREP_DEBUG_INFO
            data_ = nullptr;
#endif
        }

        XPlane(const XPlaneBase& base) { setBase(base); }
        void setBase(const XPlaneBase& base);
        void set_positive_from_id(int i) { assert(i >= 0); id_ = i + 1; }
        const Real* data() const { return base().data(); }
        const XPlaneBase& base() const;

		void inverse() { id_ = -id_; }
		XPlane opposite() const { return XPlane(-id_); }
        Real signd() const { return is_inverse() ? -1.0 : 1.0; }
        const Real* normal() const { return data(); }
        uint32_t absid() const { return is_valid() ? std::abs(id_) - 1 : INVALID_UINT32; }

		bool is_valid() const { return id_ != 0; }
        bool is_inverse() const { return id_ < 0; }

		// predicates
		Oriented_side orientation(const PlanePoint&) const;
        Oriented_side orientation(const cyPointT&) const;

        bool has_on(const cyPointT& p) const { return orientation(p) == ON_ORIENTED_BOUNDARY; }
        bool has_on(const PlanePoint& p) const { return orientation(p) == ON_ORIENTED_BOUNDARY; }
        bool id_equals(const XPlane& p) const { return std::abs(id_) == std::abs(p.id_); }

        bool coplanar(const XPlane& p) const { return base().coplanar(p.base()); }

#ifdef PREP_DEBUG_INFO
        void debug() { data_ = data(); }
#endif

    protected:
        XPlane(int i) : id_(i) {
#ifdef PREP_DEBUG_INFO
            debug();
#endif  
        }
        XPlane(int i, bool inv) : id_(inv?i:-i) {
#ifdef PREP_DEBUG_INFO
            debug();
#endif
        }

		int id_; // id = (realId+1) * sign

#ifdef PREP_DEBUG_INFO
        const Real* data_;
#endif
	};


	class PlaneLine
	{
	public:
		PlaneLine() {}
		PlaneLine(const XPlane& a, const XPlane& b) :
			planes_{ a, b } {
#ifdef PREP_DEBUG_INFO
            vec3_mul_cross(normal, a.data(), b.data());
            vec3_norm(normal, normal);
#endif
        }

		// @return 1 a->b, 0 a=b, -1 b->a
        int linearOrder(const XPlane& a, const XPlane& b) const;
		int linearOrderNoCheck(const XPlane& a, const XPlane& b) const;
        int linearOrder(const PlanePoint& a, const PlanePoint& b) const;
        int linearOrder(const cyPointT& a, const cyPointT& b) const;
        void makePositive(XPlane& input) const;
        Real dot(const XPlane&) const;
        XPlane pickPositiveVertical(const cyPointT& p) const; // 从plane triples中找到一个不平行于Line且和Line方向相同的面
        XPlane pickPositiveVertical(const PlanePoint& p) const; // 从plane triples中找到一个不平行于Line且和Line方向相同的面
        cyPointT approxNormal() const;

#ifdef PREP_DEBUG_INFO
    protected:
        vec3 normal;
#endif
	protected:
		XPlane planes_[2];
	};


    template <class K>
    CGAL::Point_3<K> convertToPoint(const XPlane& a, const XPlane& b, const XPlane& c)
    {
        typedef CGAL::Plane_3<K> CGALPLane;
        typedef CGAL::Point_3<K> point;
        CGALPLane p[3] = { { a.data()[0], a.data()[1], a.data()[2], a.data()[3] },
        { b.data()[0], b.data()[1], b.data()[2], b.data()[3] },
        { c.data()[0], c.data()[1], c.data()[2], c.data()[3] } };
        auto result = CGAL::intersection(p[0], p[1], p[2]);
        const point* res = boost::get<point>(&*result);
        return *res;
    }


	class PlanePoint
	{
    public:
        PlanePoint(const XPlane& a, const XPlane& b, const XPlane& c) :
            m_planes{ a, b, c } {
            //m_pos = convertToPoint<Depeck>(a, b, c);
#ifdef PREP_DEBUG_INFO
            auto res = convertToPoint<Depick>(a, b, c);
            coord[0] = res.x();
            coord[1] = res.y();
            coord[2] = res.z();
#endif
        }
	public:
        XPlane& plane(int i) { return m_planes[i]; }
        const XPlane& plane(int i) const { return m_planes[i]; }
        bool operator==(const PlanePoint& p) const;
        bool operator==(const cyPointT& p) const;

        cyPointT toVertexBased() const 
        { 
            auto tmp = convertToPoint<Depick>(m_planes[0], m_planes[1], m_planes[2]);
            return cyPointT(tmp.x(), tmp.y(), tmp.z());
        }

#ifdef PREP_DEBUG_INFO
    protected:
        vec3 coord;
#endif

	protected:
        //CGAL::Point_3<Depeck> m_pos;
		XPlane m_planes[3];
	};

    Real sign(const XPlane& p, const XPlane& q, const XPlane& input);
}