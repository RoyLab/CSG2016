#pragma once
#include "global.h"

namespace Boolean
{
	class XPlaneBase
	{
	public:
		XPlaneBase();
		XPlaneBase(Real *);
		XPlaneBase(Real, Real, Real, Real);

		template <class PointT>
		XPlaneBase(const PointT &p, const PointT&, const PointT& r);

		const Real& a() const { return m_data[0]; }
		const Real& b() const { return m_data[1]; }
		const Real& c() const { return m_data[2]; }
		const Real& d() const { return m_data[3]; }

		const Real* data() const { return m_data; }
	protected:
		Real m_data[4];
	};

	class XPlane
	{
	public:
		XPlane() : id(0) {}
		XPlane(int i) : id(i) {}
		XPlane(int i, bool inv) : id(inv?i:-i) {}

		const XPlaneBase& base() const;
		bool isInverse() const { return id < 0; }
		void inverse() { id = -id; }
		XPlane opposite() const { return XPlane(-id); }

		const Real* data() const { return base().data(); }
	protected:
		int id; // id = (realId+1) * sign
	};

	class XPoint
	{
	public:
		XPoint();

	protected:
		XPlane m_planes[3];
	};
}