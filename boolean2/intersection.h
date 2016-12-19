#pragma once
#include <vector>
#include <list>
#include "global.h"
#include "preps.h"
#include "Octree.h"
#include "xmemory.h"

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

    struct PBIRep
    {
		XPlane pends[2];
		uint32_t ends[2];

		std::vector<NeighborInfo> neighbor;

        virtual int is_derived_cls() const { return false; }
		virtual ~PBIRep() {}
    };

	typedef PBIRep EdgePBI;

	struct FacePBI: public PBIRep
	{
        int is_derived_cls() const { return true; }
        XPlane vertPlane;
	};

    class EdgeInsctData
    {
    public:
        typedef std::list<EdgePBI> PBIList;
        typedef std::map<uint32_t, PBIList> PBILists;
		typedef std::vector<VertexIndex> VertexList;

		void refine(void* pData);
        bool isRefined() const { return bRefined; }
        uint32_t* point(const PlanePoint&);

	public:
        class PbiPairIterator
        {
        public:
            PbiPairIterator(EdgeInsctData* target);
            void operator++();
            operator bool() const;

        private:
            PBILists::iterator outer_end_, inner_end_;
            PBILists::iterator outer_cur_, inner_cur_;

            PBIList::iterator outer_pbi_end_, inner_pbi_end_;
            PBIList::iterator outer_pbi_cur_, inner_pbi_cur_;
        };

        class PbiIterator
        {
        public:
            PbiIterator(EdgeInsctData* target);
            void operator++();
            EdgePBI* operator->() const;
            EdgePBI* pointer() const;
            operator bool() const;

        private:
            PBILists::iterator end_;
            PBILists::iterator cur_;

            PBIList::iterator pbi_end_;
            PBIList::iterator pbi_cur_;
        };

        PbiPairIterator pbi_pair_begin()
        {
            return PbiPairIterator(this);
        }

        PbiIterator pbi_begin()
        {
            return PbiIterator(this);
        }


        PBILists	inscts;
		VertexList  points;

    protected:
        bool		bRefined = false;
    };

    class FaceInsctData
    {
    public:
        // eId: -1 means tess intersection, -2 means from by propagate at that stage
        struct Vertex { VertexIndex vId; EdgeSIndex eId; };
        typedef std::list<FacePBI> PBIList;
        typedef std::map<uint32_t, PBIList> PBILists;
        typedef std::vector<Vertex> VertexList;

        void refine(void* pData);
        bool isRefined() const { return bRefined; }
        uint32_t* point(const PlanePoint&, EdgeSIndex eIdx);

    public:
        class PbiPairIterator
        {
        public:
            PbiPairIterator(FaceInsctData* target);
            void operator++();
            operator bool() const;

        private:
            PBILists::iterator outer_end_, inner_end_;
            PBILists::iterator outer_cur_, inner_cur_;

            PBIList::iterator outer_pbi_end_, inner_pbi_end_;
            PBIList::iterator outer_pbi_cur_, inner_pbi_cur_;
        };

        class PbiIterator
        {
        public:
            PbiIterator(FaceInsctData* target);
            void operator++();
            FacePBI* operator->() const;
            EdgePBI* pointer() const;
            operator bool() const;

        private:
            PBILists::iterator end_;
            PBILists::iterator cur_;

            PBIList::iterator pbi_end_;
            PBIList::iterator pbi_cur_;
        };

        PbiPairIterator pbi_pair_begin()
        {
            return PbiPairIterator(this);
        }

        PbiIterator pbi_begin()
        {
            return PbiIterator(this);
        }

        PBILists    inscts;
        VertexList  points;

    protected:
        void resolveIntersection(Triangle* pTri, std::vector<VertexIndex>* strayVertices = nullptr);
        bool		bRefined = false;
    };
}