#include "precompile.h"

#include "hybrid_geometry.h"

namespace Boolean
{
    int linear_order(const PlaneLine& l, const MyVertex& a, const MyVertex& b)
    {
        int type = 0;
        if (a.isPlaneRep()) type += 1;
        if (b.isPlaneRep()) type += 2;

        switch (type)
        {
        case 0:
            return l.linear_order(a.point(), b.point());
        case 1:
            return l.linear_order(a.ppoint(), b.point());
        case 2:
            return l.linear_order(a.point(), b.ppoint());
        case 3:
            return l.linear_order(a.ppoint(), b.ppoint());
        default:
            throw std::exception();
        }
        return 0;
    }

    XPlane pick_positive_vertical_plane(const PlaneLine &l, const MyVertex & v)
    {
        if (v.isPlaneRep())
            return l.pick_positive_vertical_plane(v.ppoint());
        else
            return l.pick_positive_vertical_plane(v.point());
    }

}