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

		virtual ~PBIRep() {}
    };

	typedef PBIRep EdgePBI;

	struct FacePBI: public PBIRep
	{
		XPlane vertPlane;
        //std::list<Triangle*> pTris; // ≤ªø…÷ÿ∏¥
	};

    class EdgeInsctData
    {
    public:
		typedef std::map<uint32_t, std::list<EdgePBI>> PBIList;
		typedef std::vector<MyVertex::Index> VertexList;

		void refine(void* pData);
        bool isRefined() const { return bRefined; }
        uint32_t* point(const XPoint&);

	public:
		PBIList		inscts;
		VertexList  points;

    protected:
        bool		bRefined = false;
    };

    class FaceInsctData
    {
    public:
        struct Vertex { MyVertex::Index vId; MyEdge::SIndex eId; };
        typedef std::map<uint32_t, std::list<FacePBI>> PBIList;
        typedef std::vector<Vertex> VertexList;

        void refine(void* pData);
        bool isRefined() const { return bRefined; }
        uint32_t* point(const XPoint&, MyEdge::SIndex eIdx);

    public:
        PBIList		inscts;
        VertexList  points;

    protected:
        void resolveIntersection(Triangle* pTri);
        bool		bRefined = false;
    };

    void doIntersection(std::vector<RegularMesh*>&, std::vector<Octree::Node*>&);
}