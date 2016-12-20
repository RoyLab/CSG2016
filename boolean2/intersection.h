#pragma once
#include <vector>

#include "preps.h"
#include "xmemory.h"
#include "hybrid_geometry.h"

namespace Boolean
{
    struct NeighborInfo
	{
		enum Type {Vertex, Edge, Face};

		Type type;
		uint32_t neighborMeshId;
		union
		{
			int neighborEdgeId; // >0 is fid, < 0 is eid,  0 is invalid
			Triangle* pTrangle;
		};
	};

    struct PbiRep
    {
		uint32_t ends[2];
		std::vector<NeighborInfo> neighbor;

        virtual int is_derived_cls() const { return false; }
		virtual ~PbiRep() {}
    };

	typedef PbiRep EdgePbi;

	struct FacePbi: public PbiRep
	{
        int is_derived_cls() const { return true; }
        XPlane vertPlane;
        XPlane pends[2];
    };

    class EdgeInsctData
    {
    public:
        typedef std::vector<EdgePbi> PbiList;
        typedef std::map<MeshIndex, PbiList> PbiLists;

        struct Vertex
        {
            VertexIndex vertex_idx;
            XPlane      plane_rep;
        };
		typedef std::vector<Vertex> VertexList;

		void refine(Triangle* triangle, int which_edge);
        bool isRefined() const { return bRefined; }
        VertexIndex* point(const PlanePoint&, const XPlane& plane);

	public:
        EdgeInsctData(const XPlane& sp, const XPlane& bp)
        {
            line = PlaneLine(sp, bp);
            line.inverse();
        }
        class PbiPairIterator
        {
        public:
            PbiPairIterator(EdgeInsctData* target);
            void operator++();
            operator bool() const;

        private:
            PbiLists::iterator outer_end_, inner_end_;
            PbiLists::iterator outer_cur_, inner_cur_;

            PbiList::iterator outer_pbi_end_, inner_pbi_end_;
            PbiList::iterator outer_pbi_cur_, inner_pbi_cur_;
        };

        class PbiIterator
        {
        public:
            PbiIterator(EdgeInsctData* target);
            void operator++();
            EdgePbi* operator->() const;
            EdgePbi* pointer() const;
            operator bool() const;

        private:
            PbiLists::iterator end_;
            PbiLists::iterator cur_;

            PbiList::iterator pbi_end_;
            PbiList::iterator pbi_cur_;
        };

        PbiPairIterator pbi_pair_begin()
        {
            return PbiPairIterator(this);
        }

        PbiIterator pbi_begin()
        {
            return PbiIterator(this);
        }

        PbiLists	inscts;
		VertexList  points;
        PlaneLine   line;

    protected:
        bool		bRefined = false;
    };

    class FaceInsctData
    {
    public:
        // eId: -1 means tess intersection, -2 means from by propagate at that stage
        struct Vertex { VertexIndex vId; EdgeSIndex eId; };
        typedef std::vector<FacePbi> PbiList;
        typedef std::map<MeshIndex, PbiList> PbiLists;
        typedef std::vector<Vertex> VertexList;

        void refine(Triangle*);
        bool isRefined() const { return bRefined; }
        VertexIndex* point(const PlanePoint&, EdgeSIndex eIdx);

    public:
        class PbiPairIterator
        {
        public:
            PbiPairIterator(FaceInsctData* target);
            void operator++();
            operator bool() const;

        private:
            PbiLists::iterator outer_end_, inner_end_;
            PbiLists::iterator outer_cur_, inner_cur_;

            PbiList::iterator outer_pbi_end_, inner_pbi_end_;
            PbiList::iterator outer_pbi_cur_, inner_pbi_cur_;
        };

        class PbiIterator
        {
        public:
            PbiIterator(FaceInsctData* target);
            void operator++();
            FacePbi* operator->() const;
            FacePbi* pointer() const;
            operator bool() const;

        private:
            PbiLists::iterator end_;
            PbiLists::iterator cur_;

            PbiList::iterator pbi_end_;
            PbiList::iterator pbi_cur_;
        };

        PbiPairIterator pbi_pair_begin()
        {
            return PbiPairIterator(this);
        }

        PbiIterator pbi_begin()
        {
            return PbiIterator(this);
        }

        PbiLists    inscts;
        VertexList  points;

    protected:
        void removeOverlapPbi();
        void resolveIntersection(Triangle*);
        bool		bRefined = false;
    };
}