#pragma once
#include <cstdint>
#include <map>

#include <macros.h>

#include "csgdefs.h"
#include "preps.h"

/// this file is strongly coupled with xmemory.h 
namespace Boolean
{
    typedef uint32_t Index;
    //typedef Index VertexIndex;
    struct VertexIndex
    {
        int b;
    };
    typedef Index EdgeIndex;
    typedef Index MeshIndex;

    typedef int32_t SIndex;
    typedef Index VertexSIndex;
    typedef Index EdgeSIndex;
    typedef Index MeshSIndex;

    class MyVertex
    {
        friend class GlobalData;
    public:
        bool isPlaneRep() const { return rep_ < 0; }
        bool isValid() const { return rep_ != 0; }

        bool isCoincident(const PlanePoint& p) const;
        bool isCoincident(const cyPointT& p) const;
        bool id_equals(const MyVertex& p) const;
        bool has_on(const XPlane& p) const;
        bool has_on(const PlaneLine& p) const
        {
            return has_on(p.plane(0)) &&
                has_on(p.plane(1));
        }

        Oriented_side orientation(const XPlane& plane) const;

        bool findEdge(EdgeIndex other, EdgeIndex* result = nullptr) const;
        const std::vector<EdgeIndex>& edges() const;
        void add_edge(EdgeIndex edge_idx);

        const PlanePoint& plane_rep() const;
        const cyPointT& vertex_rep() const;

    private:
        // only for global object usage
        void setAsPRep(int i);
        void setAsVRep(int i);

        VertexIndex absId() const { assert(rep_);  return std::abs(rep_) - 1; }

        // + is v-base, - is p-base
        int rep_ = 0;
        int merge_ = -1;
        union
        {
            ExternPtr const PlanePoint* ppoint_;
            ExternPtr const cyPointT* point_;
        };
        std::vector<EdgeIndex> edges_;
    };

    struct EdgeFaceHandle
    {
        int orientation; // +1 same, -1 oppo
        IPolygon* ptr;
    };

    class MyEdge
    {
    public:
        class ConstFaceIterator
        {
        public:
            ConstFaceIterator(const MyEdge& edge) :
                m_edge(edge), stage(0) {}

            ConstFaceIterator& operator++();
            operator bool() const { return stage != -1; }

            const EdgeFaceHandle& faceHandle() const
            {
                if (stage < 0) throw std::exception();
                if (stage < 2) return m_edge.fhs[stage];
                else return *eItr;
            }
            const IPolygon* face() const { return faceHandle().ptr; }
            int orientation() const { return faceHandle().orientation; }

        private:
            const MyEdge& m_edge;
            int stage; // > 2, then eItr
            std::vector<EdgeFaceHandle>::const_iterator eItr;
        };

        class FaceIterator
        {
        public:
            FaceIterator(MyEdge& edge, bool triangle = false);

            FaceIterator& operator++();
            FaceIterator& incrementToTriangle();

            operator bool() const { return stage != -1; }
            EdgeFaceHandle& faceHandle() const
            {
                if (stage < 0) throw std::exception();
                if (stage < 2) return m_edge.fhs[stage];
                else return *eItr;
            }
            IPolygon* face() const { return faceHandle().ptr; }
            Triangle* as_triangle() const;
            int orientation() const { return faceHandle().orientation; }

        private:
            MyEdge& m_edge;
            int stage; // > 2, then eItr
            std::vector<EdgeFaceHandle>::iterator eItr;
        };

    public:
        VertexIndex ends[2];
        std::map<MeshIndex, NeighborInfo>* neighbor = nullptr;
        EdgeInsctData* inscts = nullptr;
        //bool noOverlapNeighbor = false;

    public:
        MyEdge(VertexIndex a, VertexIndex b) : ends{ a, b }
        {
            fhs[0].ptr = 0;
            fhs[1].ptr = 0;
        }

        ~MyEdge();
        void addAjacentFace(VertexIndex s, VertexIndex e, IPolygon* fPtr);
        int faceOrientation(const IPolygon*) const;
        bool remove(const IPolygon*); // BUG: vertex √ª¥¶¿Ì
        uint32_t faceCount() const;
        VertexIndex theOtherVId(VertexIndex thiz) const;
        MyVertex& theOtherVertex(VertexIndex thiz) const;

    private:
        EdgeFaceHandle fhs[2];
        std::vector<EdgeFaceHandle> extrafhs;
    };
}