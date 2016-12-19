#pragma once
#include "preps.h"
#include "xmemory.h"

namespace Boolean
{
    // +: a->b, 0: a==b, -: a<-b
    int linear_order(const PlaneLine& l, const MyVertex& a, const MyVertex& b);
    Oriented_side orientation(const XPlane& p, const MyVertex& v);
    XPlane pickPositiveVertical(const PlaneLine&, const MyVertex&);


    class LinOrderObj
    {
    public:
        struct Item
        {
            const MyVertex* vertex;
            XPlane plane;
            int id;
        };

        bool operator()(const Item& a, const Item& b) const
        {
            int type = 0;
            if (a.vertex->isPlaneRep()) type += 1;
            if (b.vertex->isPlaneRep()) type += 2;

            switch (type)
            {
            case 0:
                return line_.linear_order(a.vertex->point(),
                    b.vertex->point()) > 0;
            case 1:
                return a.plane.orientation(b.vertex->point()) == ON_POSITIVE_SIDE;
            case 2:
                return b.plane.orientation(a.vertex->point()) == ON_NEGATIVE_SIDE;
            case 3:
                return line_.linear_order(a.plane, b.plane) > 0;
            default:
                throw std::exception();
            }
        }

    private:
        PlaneLine line_;
    };


    /// clockwise, to correctly sort, we have to use quicksort
    struct CircularOrderObj
    {
    public:
        struct Item
        {
            XPlane prep;
            int id;
        };

        CircularOrderObj(const XPlane& p) : supporting_plane_(p) {}

        bool operator() (const Item& i, const Item& j) const
        {
            Real res = sign(supporting_plane_, i.prep, j.prep);
            assert(res != Real(0));
            return res < Real(0);
        }

        bool checkSeq(std::vector<Item>& vec) const
        {
            for (int i = 1; i < vec.size(); i++)
            {
                if (!(*this)(vec[i - 1], vec[i]) && (*this)(vec[i], vec[i - 1]))
                    return false;
            }
            return true;
        }

    private:
        XPlane supporting_plane_;
    };

    enum LineInsctResultType
    {
        No_Intersect = 0, Colinear, Intersect
    };

    //struct LineInsctResult
    //{
    //    VertexType data_type;
    //    PlanePoint ppoint;
    //    cyPointT   point;
    //};

    template <class T1, class T2, class T3, class T4>
    class IResolveLineIntersection
    {
    public:
        virtual void resolve(Oriented_side side[4], const XPlane& a, const XPlane&b,
            const T1& a0, const T2& a1, const T3& b0, const T4& b1) = 0;
    };

    template <class T1, class T2, class T3, class T4>
    LineInsctResultType plane_based_line_intersection(const XPlane& a, const XPlane&b,
        const T1& a0, const T2& a1, const T3& b0, const T4& b1,
        IResolveLineIntersection<T1, T2, T3, T4>* resolve)
    {
        Oriented_side side[4] = { 0 };
        side[0] = orientation(a, b0);
        side[1] = orientation(a, b1);
        if (side[0] == side[1])
        {
            if (side[0] == ON_ORIENTED_BOUNDARY)
            {
                return Colinear;
            }
            else
            {
                return No_Intersect;
            }
        }

        side[2] = orientation(b, a0);
        side[3] = orientation(b, a1);
        if (side[2] == side[3])
        {
            assert(side[2] != ON_ORIENTED_BOUNDARY);
            return No_Intersect;
        }
        resolve->resolve(side, a, b, a0, a1, b0, b1);
        return Intersect;
    }

}