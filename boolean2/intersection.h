#pragma once
#include <vector>
#include <list>
#include "global.h"
#include "preps.h"
#include "Octree.h"
#include "xmemory.h"
#include "RegularMesh.h"

namespace Boolean
{
    struct NeighborInfo
	{
		enum Type {Vertex, Edge, Face};

		Type type;
		RegularMesh::Index neighborMeshId;
		union
		{
			int neighborEdgeId; // >0 is fid, < 0 is eid,  0 is invalid
			Triangle* pTrangle;
		};
	};

    struct PBIBase
    {
        XPlane pends[2];
        Triangle::LocalVertexId ends[2];
        std::vector<NeighborInfo> neighbor;
    };

    class EdgeInsctData
    {
    public:
        typedef PBIBase PBI;
        typedef std::map<RegularMesh::Index, std::list<PBI>> PBIList;
        typedef std::vector<MyVertex::Index> VertexList;

		void refine(void* pData);
        bool isRefined() const { return bRefined; }
        uint32_t* point(const XPoint&);
        Triangle::LocalVertexId localId(const XPoint&, MyVertex::Index *&slot);
	public:
		PBIList		inscts;
		VertexList  points;

    protected:
        bool		bRefined = false;
    };

    class FaceInsctData
    {
    public:
        struct PBI: public PBIBase
        {
            XPlane vertPlane;
        };

        struct Vertex { MyVertex::Index vId; MyEdge::SIndex eId; };
        typedef std::map<RegularMesh::Index, std::list<PBI>> PBIList;
        typedef std::vector<Vertex> VertexList;

        void refine(void* pData);
        bool isRefined() const { return bRefined; }
        uint32_t* point(const XPoint&, MyEdge::SIndex eIdx);
        Triangle::LocalVertexId localId(const XPoint&, MyEdge::SIndex eIdx, MyVertex::Index *&slot);

    public:
        PBIList		inscts;
        VertexList  points;
        std::set<RegularMesh::Index> insctMeshes;

    protected:
        void resolveIntersection(Triangle* pTri, std::vector<MyVertex::Index>* strayVertices = nullptr);
        bool bRefined = false;
    };

    void doIntersection(std::vector<RegularMesh*>&, std::vector<Octree::Node*>&);
}