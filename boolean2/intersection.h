#pragma once
#include <vector>

#include <xstruct.hpp>

#include "preps.h"
#include "xmemory.h"
#include "hybrid_geometry.h"

namespace Boolean
{
    struct NeighborInfo
	{
		enum Type {Vertex, Edge, Face};

		Type type;
		union
		{
			int neighborEdgeId;
			Triangle* pTrangle;
		};
	};

    struct PbiRep
    {
		uint32_t ends[2];
		std::map<MeshIndex, NeighborInfo> neighbor;

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
        VertexIndex* find_point(const PlanePoint & p);
        VertexIndex* find_point(const XPlane & p);

        // debug
        bool checkOrientation(const EdgePbi* pbi) const
        {
            return linear_order(line, pbi->ends[0], pbi->ends[1]) > 0;
        }

        bool checkOrientation() const
        {
            //check pbi
            auto itr = const_cast<EdgeInsctData*>(this)->pbi_begin(); // lazy to write a const version of iterator
            for (; itr; ++itr)
            {
                if (!checkOrientation(itr.pointer()))
                {
                    return false;
                }
            }
            // check point
            for (auto& p : points)
            {
                if (line.dot(p.plane_rep) <= 0)
                {
                    return false;
                }

                if (!xvertex(p.vertex_idx).has_on(line))
                {
                    return false;
                }
            }
            return true;
        }

	public:
        EdgeInsctData(const XPlane& sp, const XPlane& bp, const IPolygon* poly, const MyEdge* edge)
        {
            line = PlaneLine(sp, bp);
            int face_ori = edge->faceOrientation(poly);
            assert(face_ori != 0);
            if (face_ori > 0)
            {
                line.inverse();
            }
            assert(linear_order(line, edge->ends[0], edge->ends[1]) > 0);
        }
        //class PbiPairIterator
        //{
        //public:
        //    PbiPairIterator(EdgeInsctData* target);
        //    void operator++();
        //    operator bool() const;

        //private:
        //    PbiLists::iterator outer_end_, inner_end_;
        //    PbiLists::iterator outer_cur_, inner_cur_;

        //    PbiList::iterator outer_pbi_end_, inner_pbi_end_;
        //    PbiList::iterator outer_pbi_cur_, inner_pbi_cur_;
        //};

        //class PbiIterator
        //{
        //public:
        //    PbiIterator(EdgeInsctData* target);
        //    void operator++();
        //    EdgePbi* operator->() const;
        //    EdgePbi* pointer() const;
        //    operator bool() const;

        //private:
        //    void assignInnerPtr();

        //    PbiLists::iterator end_;
        //    PbiLists::iterator cur_;

        //    PbiList::iterator pbi_end_;
        //    PbiList::iterator pbi_cur_;
        //};

        //PbiPairIterator pbi_pair_begin()
        //{
        //    return PbiPairIterator(this);
        //}

        static void assignInnerPtr(PbiLists::iterator o, 
            PbiList::iterator& begin, PbiList::iterator& end)
        {
            begin = o->second.begin();
            end = o->second.end();
        }

        typedef XR::DoubleIterator<
            PbiLists::iterator, 
            PbiList::iterator,
            assignInnerPtr>
            PbiIterator;

        PbiIterator pbi_begin()
        {
            return PbiIterator(inscts.begin(), inscts.end());
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
        struct Vertex 
        { 
            VertexIndex vId; 
            EdgeSIndex eId;
            PlaneLine line;
        };
        typedef std::vector<FacePbi> PbiList;
        typedef std::map<MeshIndex, PbiList> PbiLists;
        typedef std::vector<Vertex> VertexList;

        FaceInsctData(const Triangle* triangle) : splane_(triangle->supportingPlane()) {}

        void refine(Triangle*);
        bool isRefined() const { return bRefined; }
        Vertex* point(const PlanePoint&, EdgeSIndex eIdx, const PlaneLine&);

        bool checkOrientation(const FacePbi* pbi) const;
        bool checkOrientation(const Vertex& vertex) const;
        bool checkOrientation(const PlanePoint& vertex, EdgeSIndex eId) const;
        bool checkOrientation(const PlanePoint& vertex) const;
        bool checkOrientation() const
        {
            //check pbi
            auto itr = const_cast<FaceInsctData*>(this)->pbi_begin(); // lazy to write a const version of iterator
            for (; itr; ++itr)
            {
                if (!checkOrientation(itr.pointer()))
                {
                    return false;
                }
            }
            // check point
            for (auto& p : points)
            {
                if (!checkOrientation(p))
                {
                    return false;
                }
            }
            return true;
        }

    public:

        static void assignInnerPtr(PbiLists::iterator o,
            PbiList::iterator& begin, PbiList::iterator& end)
        {
            begin = o->second.begin();
            end = o->second.end();
        }

        typedef XR::DoubleIterator<
            PbiLists::iterator,
            PbiList::iterator,
            assignInnerPtr>
            PbiIterator;

        PbiIterator pbi_begin()
        {
            return PbiIterator(inscts.begin(), inscts.end());
        }

        PbiLists    inscts;
        VertexList  points;

    protected:
        void removeOverlapPbi();
        void resolveIntersection(Triangle*);
        bool		bRefined = false;
        XPlane splane_;
    };
}