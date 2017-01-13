#pragma once
#ifdef PREP_DEBUG_INFO
#include <CGAL/intersection_3.h>
#endif
#include <macros.h>
#include "global.h"

namespace Boolean
{
    typedef Real Real4[4];

    class XPlaneBase
    {
        friend class XPlane;
    public:
        XPlaneBase() {}
        XPlaneBase(Real x, Real y, Real z, Real w) :
            data_{ x,y,z,w } {}
        //virtual ~XPlaneBase() {}

        const Real& a() const { return data_[0]; }
        const Real& b() const { return data_[1]; }
        const Real& c() const { return data_[2]; }
        const Real& d() const { return data_[3]; }
        const Real* data() const { return data_; }
        const Real& operator[](int i) const { return data_[i]; }
        const Real& at(int i) const { return data_[i]; }

    private:
        Real& operator[](int i) { return data_[i]; }
        Real& at(int i) { return data_[i]; }

        Real data_[4];
    };

	class XPlane
	{
	public:
		XPlane() : id_(0)
#ifdef PREP_DEBUG_INFO
            ,debug_data_(nullptr)
#endif
        {}

        //XPlane(const XPlaneBase& base) { setBase(base); }

        void construct_from_three_vertices(const cyPointT& p, const cyPointT& q, const cyPointT& r);
        void construct_from_one_vertex_two_edges(const cyPointT& p, const cyPointT& q, const cyPointT& r);
        void construct_coicident_plane(const PlaneLine & l, const cyPointT & p);

        // state
        bool is_valid() const { return id_ != 0; }
        bool is_inverse() const { return id_ < 0; }
        bool is_degenerate() const;

        // access
        const XPlaneBase& get_base() const;
        const Real* get_data() const { return get_base().data(); }
        const Real* normal() const { return get_data(); }
        uint32_t absid() const { return is_valid() ? std::abs(id_) - 1 : INVALID_UINT32; }
        Real signd() const { return is_inverse() ? -1.0 : 1.0; }

        //void set_positive_from_id(int i) { assert(i >= 0); id_ = i + 1; }

        // manipulate
		void inverse() { id_ = -id_; }
        XPlane opposite() const { XPlane res; res.id_ = -id_; return res; }

		// predicates
		Oriented_side orientation(const PlanePoint&) const;
        Oriented_side orientation(const cyPointT&) const;
        bool parallel(const XPlane& p) const;
        bool normal_equals(const XPlane& p) const;
        bool has_on(const cyPointT& p) const { return orientation(p) == ON_ORIENTED_BOUNDARY; }
        bool has_on(const PlanePoint& p) const { return orientation(p) == ON_ORIENTED_BOUNDARY; }
        bool id_equals(const XPlane& p) const { return std::abs(id_) == std::abs(p.id_); }

    protected:
//        XPlane(int i) : id_(i) {
//#ifdef PREP_DEBUG_INFO
//            debug();
//#endif  
//        }
//        XPlane(int i, bool inv) : id_(inv?i:-i) {
//#ifdef PREP_DEBUG_INFO
//            debug();
//#endif
        //}
        XPlaneBase* register_base();

		int id_; // id = (realId+1) * sign

#ifdef PREP_DEBUG_INFO
        ExternPtr const Real (*debug_data_)[4];
#endif
	};

    Real orientation(const XPlane& p, const XPlane& q, const XPlane& input);
    Real orientation(const XPlane & p, const XPlane & q, const XPlane & r, const XPlane & s);

	class PlaneLine
	{
	public:
		PlaneLine() {}
		PlaneLine(const XPlane& a, const XPlane& b) :
			planes_{ a, b } {
            assert(!a.parallel(b));
#ifdef PREP_DEBUG_INFO
            //vec3_mul_cross(normal, a.get_data(), b.get_data());
            //vec3_norm(normal, normal);
            cyPointT approx = approxNormal();
            vec3_copy(normal, reinterpret_cast<Real*>(&approx));
#endif
        }

		// @return 1 a->b, 0 a=b, -1 b->a
        int linear_order(const XPlane& a, const XPlane& b) const;
        int linear_order(const PlanePoint& a, const XPlane& b) const;
        int linear_order(const XPlane& a, const PlanePoint& b) const;

        int linear_order(const cyPointT& a, const cyPointT& b) const;
        int linear_order(const PlanePoint& a, const cyPointT& b) const;
        int linear_order(const cyPointT& a, const PlanePoint& b) const;

        int linear_order(const XPlane& a, const cyPointT& b) const;
        int linear_order(const cyPointT& a, const XPlane& b) const;

        int linear_order(const PlanePoint& a, const PlanePoint& b) const;

        int linear_order_unsafe(const XPlane& a, const XPlane& b) const;
        int linear_order_unsafe(const XPlane& a, const PlanePoint& b) const;
        int linear_order_unsafe(const PlanePoint& a, const XPlane& b) const;

        int linear_order_unsafe(const XPlane& a, const cyPointT& b) const;
        int linear_order_unsafe(const cyPointT& a, const XPlane& b) const;

        bool linear_coincident_no_check(const XPlane& a, const XPlane& b) const;

        void make_positive(XPlane& input) const;
        Real dot(const XPlane&) const;
        Real dot(const Real*) const;
        XPlane pick_positive_vertical_plane(const cyPointT& p) const; // 从plane triples中找到一个不平行于Line且和Line方向相同的面
        XPlane pick_positive_vertical_plane(const PlanePoint& p) const; // 从plane triples中找到一个不平行于Line且和Line方向相同的面
        void inverse();
        cyPointT approxNormal() const;
        const XPlane& plane(int i) const { return planes_[i]; }

        template <class Point>
        bool has_on(const Point& p) const
        {
            return planes_[0].has_on(p) && 
                planes_[1].has_on(p); 
        }

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
        CGALPLane p[3] = { { a.get_data()[0], a.get_data()[1], a.get_data()[2], a.get_data()[3] },
        { b.get_data()[0], b.get_data()[1], b.get_data()[2], b.get_data()[3] },
        { c.get_data()[0], c.get_data()[1], c.get_data()[2], c.get_data()[3] } };
        auto result = CGAL::intersection(p[0], p[1], p[2]);
        const point* res = boost::get<point>(&*result);
        return *res;
    }


	class PlanePoint
	{
    public:
        PlanePoint(const XPlane& a, const XPlane& b, const XPlane& c);
	public:
        XPlane& plane(int i) { return planes_[i]; }
        const XPlane& plane(int i) const { return planes_[i]; }
        bool value_equals(const PlanePoint& p) const;
        bool value_equals(const cyPointT& p) const;

        bool check_positive() const;

        cyPointT toVertexBased() const 
        { 
            auto tmp = convertToPoint<Depick>(planes_[0], planes_[1], planes_[2]);
            return cyPointT(tmp.x(), tmp.y(), tmp.z());
        }

#ifdef PREP_DEBUG_INFO
    protected:
        vec3 coord;
#endif

	protected:
        //CGAL::Point_3<Depeck> m_pos;
		XPlane planes_[3];
	};
}