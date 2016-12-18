#pragma once
#include "preps.h"
#include "xmemory.h"

namespace Boolean
{
    // +: a->b, 0: a==b, -: a<-b
    int linearOrder(const PlaneLine& l, const MyVertex& a, const MyVertex& b);
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
                return line_.linearOrder(a.vertex->point(),
                    b.vertex->point()) > 0;
            case 1:
                return a.plane.orientation(b.vertex->point()) == ON_POSITIVE_SIDE;
            case 2:
                return b.plane.orientation(a.vertex->point()) == ON_NEGATIVE_SIDE;
            case 3:
                return line_.linearOrder(a.plane, b.plane) > 0;
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
}