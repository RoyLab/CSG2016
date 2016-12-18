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
        //std::list<Triangle*> pTris; // �����ظ�
	};

    class EdgeInsctData
    {
    public:
		typedef std::map<uint32_t, std::list<EdgePBI>> PBIList;
		typedef std::vector<VertexIndex> VertexList;

		void refine(void* pData);
        bool isRefined() const { return bRefined; }
        uint32_t* point(const PlanePoint&);

	public:
		PBIList		inscts;
		VertexList  points;

    protected:
        bool		bRefined = false;
    };

    class FaceInsctData
    {
    public:
        // eId: -1 means tess intersection, -2 means from by propagate at that stage
        struct Vertex { VertexIndex vId; EdgeSIndex eId; };
        typedef std::map<uint32_t, std::list<FacePBI>> PBIList;
        typedef std::vector<Vertex> VertexList;

        void refine(void* pData);
        bool isRefined() const { return bRefined; }
        uint32_t* point(const PlanePoint&, EdgeSIndex eIdx);

    public:
        PBIList		inscts;
        VertexList  points;

    protected:
        void resolveIntersection(Triangle* pTri, std::vector<VertexIndex>* strayVertices = nullptr);
        bool		bRefined = false;
    };
}