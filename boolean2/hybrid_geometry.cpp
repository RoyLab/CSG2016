#include "precompile.h"

#include "hybrid_geometry.h"

namespace Boolean
{
    int linear_order(const PlaneLine& l, const MyVertex& a, const MyVertex& b)
    {
        int type = 0;
        if (a.isPlaneRep()) type += 1;
        if (b.isPlaneRep()) type += 2;

        Oriented_side side;
        switch (type)
        {
        case 0:
            return l.linear_order(a.point(), b.point());
        case 1:
            side = l.pick_positive_vertical_plane(a.ppoint())
                .orientation(b.point());
            if (side == ON_POSITIVE_SIDE) return 1;
            else if (side == ON_NEGATIVE_SIDE) return -1;
            else return 0;
        case 2:
            side = l.pick_positive_vertical_plane(b.ppoint())
                .orientation(a.point());
            if (side == ON_NEGATIVE_SIDE) return 1;
            else if (side == ON_POSITIVE_SIDE) return -1;
            else return 0;
        case 3:
            return l.linear_order(a.ppoint(), b.ppoint());
        default:
            throw std::exception();
        }
        return 0;
    }

    Oriented_side orientation(const XPlane& p, const MyVertex& v)
    {
        if (v.isPlaneRep())
            return p.orientation(v.ppoint());
        else return p.orientation(v.point());
    }

    XPlane pickPositiveVertical(const PlaneLine &l, const MyVertex & v)
    {
        if (v.isPlaneRep())
            return l.pick_positive_vertical_plane(v.ppoint());
        else
            return l.pick_positive_vertical_plane(v.point());
    }

}