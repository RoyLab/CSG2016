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

	template <class PBI>
    class InsctData
    {
    public:
		typedef std::map<uint32_t, std::list<PBI>> PBIList;
		typedef std::list<MyVertex::Index> VertexList;

		void refine(void* pData);
        bool isRefined() const { return bRefined; }
        uint32_t* point(const XPoint&);

    protected:
        bool		bRefined = false;

	public:
		PBIList		inscts;
		VertexList  points;
    };

    void doIntersection(std::vector<RegularMesh*>&, std::vector<Octree::Node*>&);

    template<class PBI>
    inline uint32_t * InsctData<PBI>::point(const XPoint &p)
    {
        for (auto itr = points.begin(); itr != points.end(); itr++)
        {
            if (xvertex(*itr) == p)
                return &*itr;
        }

        points.push_back(INVALID_UINT32);
        return &points.back();
    }
}