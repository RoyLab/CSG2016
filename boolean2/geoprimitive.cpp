#include "precompile.h"
#include "geoprimitive.h"

#include "xmemory.h"
#include "intersection.h"

namespace Boolean
{
    bool MyVertex::findEdge_local(VertexIndex other, EdgeIndex * result) const
    {
        const std::vector<EdgeIndex>* edges = nullptr;
        //if (merge_ < 0)
        //{
        //    edges = &edges_;
        //}
        //else {
        //    edges = &GlobalData::getObject()->get_merged_edges(merge_);
        //}
        edges = &edges_;
        for (EdgeIndex item : *edges)
        {
            MyEdge &e = xedge(item);
            if (vertex_id_equals_simple(e.ends[0], other) || vertex_id_equals_simple(e.ends[1], other))
            {
                if (result)
                {
                    *result = item;
                }
                return true;
            }
        }
        return false;

    }

    const std::vector<EdgeIndex>& MyVertex::edges_local() const
    {
        return edges_;
        //if (merge_ < 0)
        //{
        //    return edges_;
        //}
        //else
        //{
        //    return GlobalData::getObject()->get_merged_edges(merge_);
        //}
    }

    void MyVertex::add_edge_local(EdgeIndex edge_idx)
    {
        edges_.push_back(edge_idx);
        //if (merge_ < 0)
        //{
        //    edges_.push_back(edge_idx);
        //}
        //else
        //{
        //    GlobalData::getObject()->add_merged_edges(merge_, edge_idx);
        //}
    }

    Oriented_side MyVertex::orientation(const XPlane & p) const
    {
        if (isPlaneRep())
        {
            return p.orientation(plane_rep());
        }
        else
        {
            return p.orientation(vertex_rep());
        }
    }

    bool MyVertex::isCoincident(const PlanePoint & p) const
    {
        if (isPlaneRep())
        {
            return p.value_equals(plane_rep());
        }
        else
        {
            return p.value_equals(vertex_rep());
        }
    }

    bool MyVertex::isCoincident(const cyPointT & p) const
    {
        if (isPlaneRep())
        {
            return plane_rep().value_equals(p);
        }
        else
        {
            return (p == vertex_rep()) == 0 ? false: true;
        }
    }

    bool MyVertex::merge_equals(const MyVertex & p) const
    {
        return (merge_ > -1 && merge_ == p.merge_);
    }

    bool MyVertex::has_on(const XPlane & p) const
    {
        if (isPlaneRep())
        {
            return p.has_on(plane_rep());
        }
        else
        {
            return p.has_on(vertex_rep());
        }
    }

    const PlanePoint & MyVertex::plane_rep() const
    {
#ifdef XR_DEBUG
        if (!isValid() || !isPlaneRep())
            throw std::exception();
#endif
        //return xppoint(absId());
        return *ppoint_;
    }

    const cyPointT & MyVertex::vertex_rep() const
    {
#ifdef XR_DEBUG
        if (!isValid() || isPlaneRep())
            throw std::exception();
#endif
        //return xpoint(absId());
        return *point_;
    }

    void MyVertex::setAsPRep(int i)
    {
        rep_ = -(i + 1);
        ppoint_ = &xcppoint(i);
    }

    void MyVertex::setAsVRep(int i)
    {
        rep_ = (i + 1);
        point_ = &xcpoint(i);
    }

    ////////////////////////////////////////////////////////////////

    MyEdge::~MyEdge()
    {
        SAFE_DELETE(inscts);
        SAFE_DELETE(neighbor);
    }

    void MyEdge::addAjacentFace(VertexIndex s, VertexIndex e, IPolygon * fPtr)
    {
        assert(s == ends[0] && e == ends[1] 
            || s == ends[1] && e == ends[0]);

        int ori = 1;
        if (!vertex_id_equals_simple(s, ends[0])) ori = -1;

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

        extrafhs.emplace_back(EdgeFaceHandle{ ori, fPtr });
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
        if (vertex_id_equals_simple(ends[0], thiz))
        {
            return ends[1];
        }
        else
        {
            //assert(ends[1] == thiz);
            assert(vertex_id_equals_simple(ends[1], thiz));
            return ends[0];
        }
    }

    MyEdge::ConstFaceIterator & MyEdge::ConstFaceIterator::operator++()
    {
        if (stage == -1) throw std::exception();

        if (stage < 2)
        {
            ++stage;
            if (stage == 2)
            {
                if (m_edge.extrafhs.empty())
                {
                    stage = -1;
                    return *this;
                }
                else eItr = m_edge.extrafhs.begin();
            }
        }
        else
        {
            ++eItr;
            if (eItr == m_edge.extrafhs.end())
                stage = -1;
        }
        return *this;
    }

    MyEdge::FaceIterator::FaceIterator(MyEdge & edge, bool triangle) :
        m_edge(edge), stage(0)
    {
        if (triangle)
        {
            if (face()->getType() != IPolygon::TRIANGLE)
                incrementToTriangle();
        }
    }

    MyEdge::FaceIterator & MyEdge::FaceIterator::operator++()
    {
        if (stage == -1) throw std::exception();

        if (stage < 2)
        {
            ++stage;
            if (stage == 2)
            {
                if (m_edge.extrafhs.empty())
                {
                    stage = -1;
                    return *this;
                }
                else eItr = m_edge.extrafhs.begin();
            }
        }
        else
        {
            ++eItr;
            if (eItr == m_edge.extrafhs.end())
                stage = -1;
        }
        return *this;
    }

    MyEdge::FaceIterator & MyEdge::FaceIterator::incrementToTriangle()
    {
        do {
            ++*this;
        } while (*this && face()->getType() != IPolygon::TRIANGLE);
        return *this;
    }

    Triangle * MyEdge::FaceIterator::as_triangle() const
    {
        assert(faceHandle().ptr->getType() == IPolygon::TRIANGLE);
        return reinterpret_cast<Triangle*>(faceHandle().ptr);
    }
}