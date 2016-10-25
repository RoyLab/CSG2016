#pragma once
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
		XPlane() : id(0) {}
		XPlane(int i) : id(i) {}
		XPlane(int i, bool inv) : id(inv?i:-i) {}
        XPlane(const cyPointT& p, const cyPointT& q, const cyPointT& r);


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
        bool has_on(const cyPointT& p) { return orientation(p) == REL_ON_BOUNDARY; }

        void setId(int i) { assert(i >= 0); id = i + 1; }
        void setFromPEE(const cyPointT& p, const cyPointT& e0, const cyPointT& e1);
    protected:
		int id; // id = (realId+1) * sign
	};


	class XLine
	{
	public:
		XLine() {}
		XLine(const XPlane& a, const XPlane& b) :
			m_planes{ a, b } {}

		// @return 1 a->b, 0 a=b, -1 b->a
		int linearOrder(const XPlane& a, const XPlane& b);
		// @return 1 a->b, 0 a=b, -1 b->a
		int linearOrderNoCheck(const XPlane& a, const XPlane& b);
		void makePositive(XPlane& input);
	protected:
		XPlane m_planes[2];
	};

	class XPoint
	{
	public:
		XPoint() {}
		XPoint(const XPlane& a, const XPlane& b, const XPlane& c):
			m_planes{ a, b, c } {}

        XPlane& plane(int i) { return m_planes[i]; }
        const XPlane& plane(int i) const { return m_planes[i]; }
        bool operator==(const XPoint& p) const;
        bool operator==(const cyPointT& p) const;

	protected:
		XPlane m_planes[3];
	};



    static inline void makePositive(const XPlane& p, const XPlane& q, XPlane& input)
    {
        const Real* mat[3] = { p.data(), q.data(), input.data() };
        if (adaptiveDet3x3Sign(mat) < 0.0)
            input.inverse();
    }
}