#pragma once
#define PREP_DEBUG_INFO
#ifdef PREP_DEBUG_INFO
#include <CGAL/intersection_3.h>
#endif
#include "global.h"
#include "adaptive.h"


namespace Boolean
{
	class XPlaneBase
	{
	public:
        XPlaneBase() {}
		XPlaneBase(Real *);
		XPlaneBase(Real, Real, Real, Real);

		XPlaneBase(const cyPointT& p, const cyPointT& q, const cyPointT& r);
        void setFromPEE(const cyPointT& p, const cyPointT& e0, const cyPointT& e1);
        Oriented_side orientation(const cyPointT&) const;

		const Real& a() const { return m_data[0]; }
		const Real& b() const { return m_data[1]; }
		const Real& c() const { return m_data[2]; }
		const Real& d() const { return m_data[3]; }

        const Real* data() const { return m_data; }
        Real* data() { return m_data; }
	protected:
		Real m_data[4];
	};

	class XPlane
	{
	public:
		XPlane() : id(0) {
#ifdef PREP_DEBUG_INFO
            m_data = nullptr;
#endif
        }
		XPlane(int i) : id(i) {
#ifdef PREP_DEBUG_INFO
            debug();
#endif  
        }
		XPlane(int i, bool inv) : id(inv?i:-i) {
#ifdef PREP_DEBUG_INFO
            debug();
#endif
        }
        XPlane(const cyPointT& p, const cyPointT& q, const cyPointT& r);

#ifdef PREP_DEBUG_INFO
        void debug() { m_data = data(); }
#endif
		bool isInverse() const { return id < 0; }
		void inverse() { id = -id; }
		XPlane opposite() const { return XPlane(-id); }
        Real signd() const { return isInverse() ? -1.0 : 1.0; }
        const Real* normal() const { return data(); }

		bool isValid() const { return id != 0; }
		const XPlaneBase& base() const;
		const Real* data() const { return base().data(); }

		// predicates
		Oriented_side orientation(const XPoint&) const;
        Oriented_side orientation(const cyPointT&) const;
        bool has_on(const cyPointT& p) { return orientation(p) == ON_ORIENTED_BOUNDARY; }

        void setId(int i) { assert(i >= 0); id = i + 1; }
        void setFromPEE(const cyPointT& p, const cyPointT& e0, const cyPointT& e1);
    protected:
		int id; // id = (realId+1) * sign

#ifdef PREP_DEBUG_INFO
        const Real* m_data;
#endif
	};


	class XLine
	{
	public:
		XLine() {}
		XLine(const XPlane& a, const XPlane& b) :
			m_planes{ a, b } {
#ifdef PREP_DEBUG_INFO
            vec3_mul_cross(normal, a.data(), b.data());
            vec3_norm(normal, normal);
#endif
        }

		// @return 1 a->b, 0 a=b, -1 b->a
        int linearOrder(const XPlane& a, const XPlane& b) const;
		int linearOrderNoCheck(const XPlane& a, const XPlane& b) const;
        int linearOrder(const XPoint& a, const XPoint& b) const;
        int linearOrder(const cyPointT& a, const cyPointT& b) const;
        void makePositive(XPlane& input) const;
        Real dot(const XPlane&) const;
        XPlane pickPositiveVertical(const XPoint& p) const;

#ifdef PREP_DEBUG_INFO
    protected:
        vec3 normal;
#endif
	protected:
		XPlane m_planes[2];
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


	class XPoint
	{
	public:
		XPoint() {}
		XPoint(const XPlane& a, const XPlane& b, const XPlane& c):
			m_planes{ a, b, c } {
#ifdef PREP_DEBUG_INFO
            auto res = convertToPoint<Depick>(a, b, c);
            coord[0] = res.x();
            coord[1] = res.y();
            coord[2] = res.z();
#endif
        }

        XPlane& plane(int i) { return m_planes[i]; }
        const XPlane& plane(int i) const { return m_planes[i]; }
        bool operator==(const XPoint& p) const;
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
		XPlane m_planes[3];
	};

    Real sign(const XPlane& p, const XPlane& q, const XPlane& input);
}