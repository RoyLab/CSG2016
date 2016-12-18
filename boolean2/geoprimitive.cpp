#include "precompile.h"
#include "geoprimitive.h"
#include "xmemory.h"

namespace Boolean
{
    bool MyVertex::findEdge(uint32_t other, uint32_t * result) const
    {
        for (auto& item : edges)
        {
            auto &e = xedge(item);
            if (e.ends[0] == other || e.ends[1] == other)
            {
                if (result)
                    *result = item;
                return true;
            }
        }
        return false;
    }

    Oriented_side MyVertex::orientation(const XPlane & p) const
    {
        if (rep > 0)
        {
            return p.orientation(xpoint(rep - 1));
        }
        else
        {
            return p.orientation(xppoint((-rep) - 1));
        }
    }

    bool MyVertex::isCoincident(const PlanePoint & p) const
    {
        if (isPlaneRep())
            return xcppoint(absId()) == p;
        else
            return p == xcpoint(absId());
    }

    const PlanePoint & MyVertex::ppoint() const
    {
        if (!isValid() || !isPlaneRep())
            throw std::exception();

        return xppoint(absId());
    }

    const cyPointT & MyVertex::point() const
    {
        if (!isValid() || isPlaneRep())
            throw std::exception();

        return xpoint(absId());
    }

    MyEdge::~MyEdge()
    {
        //SAFE_DELETE(inscts);
        //SAFE_DELETE(neighbor);
    }

    void MyEdge::addAjacentFace(VertexIndex s, VertexIndex e, IPolygon * fPtr)
    {
        int ori = 1;
        if (s != ends[0]) ori = -1;

        FaceIterator itr(*this);
        for (; itr; ++itr)
        {
            if (!itr.face())
            {
                itr.faceHandle().ptr = fPtr;
                itr.faceHandle().orientation = ori;
                return;
            }
        }

        //assert(0);
        extrafhs.push_back(FH{ ori, fPtr });
    }

    int MyEdge::faceOrientation(const IPolygon * ptr) const
    {
        if (fhs[0].ptr == ptr) return fhs[0].orientation;
        if (fhs[1].ptr == ptr) return fhs[1].orientation;

        for (auto& fh : extrafhs)
            if (fh.ptr == ptr) return fh.orientation;

        return 0;
    }

    bool MyEdge::remove(const IPolygon *ptr)
    {
        FaceIterator itr(*this);
        while (itr)
        {
            if (itr.face() == ptr)
            {
                itr.faceHandle().ptr = nullptr;
                return true;
            }
            ++itr;
        }
        return false;
    }

    uint32_t MyEdge::faceCount() const
    {
        auto fItr = ConstFaceIterator(*this);
        uint32_t count = 0;
        for (; fItr; ++fItr)
        {
            if (fItr.face()->isValid())
                count++;
        }
        return count;
    }

    MyVertex & MyEdge::theOtherVertex(VertexIndex thiz) const
    {
        return xvertex(theOtherVId(thiz));
    }

    VertexIndex MyEdge::theOtherVId(VertexIndex thiz) const
    {
        if (ends[0] == thiz)
        {
            return ends[1];
        }
        else
        {
            assert(ends[1] == thiz);
            return ends[0];
        }
    }


}