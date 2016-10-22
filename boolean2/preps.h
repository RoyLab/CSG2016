#pragma once
#include "global.h"

namespace Boolean
{
    class XPlane
    {
    public:
        XPlane();
        XPlane(Real *);
        XPlane(Real x, Real y, Real z);

        template <class PointT>
        XPlane(const PointT &p, const PointT&, const PointT& r);

    protected:
        Real* m_data;
    };

    class XPoint
    {
    public:
        XPoint();
        XPoint(XPlane*);
    protected:
        XPlane *m_planes[3];
    };

}