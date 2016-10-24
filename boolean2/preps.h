#pragma once
#include "global.h"

namespace Boolean
{
    class XPlane
    {
    public:
        XPlane();
        XPlane(Real *);
        XPlane(Real, Real, Real, Real);

        template <class PointT>
        XPlane(const PointT &p, const PointT&, const PointT& r);

		const Real& a() const { return m_data[0]; }
		const Real& b() const { return m_data[1]; }
		const Real& c() const { return m_data[2]; }
		const Real& d() const { return m_data[3]; }
		
    protected:
        Real m_data[4];
    };

    class XPoint
    {
    public:
        XPoint();

    protected:
        size_t planeIds[3];
    };
	
	template <class T>

	template <class PointT, class ForwardIterator>
	void normalizeAndFilter(typename PointT::Scalar* c, typename PointT::Scalar* d,
		ForwardIterator begin, ForwardIterator end)
	{

	}
}