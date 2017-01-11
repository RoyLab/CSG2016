#pragma once
#include <macros.h>
#include <xgeometry.h>

#include <deque>

#include "global.h"
#include "preps.h"
#include "offio.h"
#include "geoprimitive.h"

namespace Boolean
{
    class IPolygon
    {
    public:
        enum TYPE {TRIANGLE, SUBPOLYGON, SUBPOLYGON_WITH_HOLES};
        int mark = UNVISITED;

    public:
		IPolygon(uint32_t i, uint32_t meshId):
			m_id(i), m_meshId(meshId) {}
		virtual ~IPolygon() {}

		uint32_t id() const { return m_id; }
        uint32_t meshId() const { return m_meshId; }

        virtual void get_vertices_for_dumping(std::vector<VertexIndex>&) const = 0;
        virtual void getAllEdges(std::vector<EdgeIndex>&) const = 0;
        virtual TYPE getType() const = 0;
        virtual VertexIndex get_rep_vertex(EdgeIndex edgeIdx) const = 0;

        //virtual MyEdge& edge(int i) const;
        //virtual uint32_t edgeId(int i) const = 0;
        //virtual MyVertex& vertex(int i)const;
        //virtual uint32_t vertexId(int i)const = 0;
        virtual bool isValid() const = 0;

        XPlane supportingPlane() const { assert(sPlane.is_valid()); return sPlane; }
        XPlane sPlane;
    protected:
        const uint32_t m_id;
        const uint32_t m_meshId;
    };

    class Triangle : public IPolygon
    {
		friend class RegularMesh;
    public:
		FaceInsctData* inscts = nullptr;
        bool add_as_insct_triangle = false;
        int octree_mark = -1;

		Triangle(uint32_t meshId, uint32_t i): IPolygon(i, meshId) {}
        ~Triangle();

		// access
		//const CGALTriangle& triangle() const { return cgalTri; }
		XPlane boundingPlane(int i) const { return bPlanes[i]; }
        cyPointT& point(int i) const;
        uint32_t degree() const { return 3; }

        MyEdge& edge(int i) const;
        MyVertex& vertex(int i)const;
        EdgeIndex edgeId(int i) const { return eIds[i]; }
        VertexIndex vertexId(int i) const { return vIds[i]; }

        void get_vertices_for_dumping(std::vector<VertexIndex>&) const;
        void getAllEdges(std::vector<EdgeIndex>& output) const;
        bool coherentEdge(int whichEdge) const; // 等价于MyEdge的faceOrientation
        VertexIndex getTheOtherVertex(EdgeIndex eId) const;
        VertexIndex get_rep_vertex(EdgeIndex edgeIdx) const;

		// search
        VertexIndex findFaceVertex(const PlanePoint& pt, EdgeIndex eIdx, const PlaneLine& line, VertexIndex*&, VertexIndex* hint = nullptr);
        VertexIndex findEdgeVertex(const PlanePoint& pt, PosTag tag, const XPlane& vert, VertexIndex*&, VertexIndex* hint = nullptr);
        VertexIndex findEdgeVertex(const PlanePoint&, PosTag, const XPlane&, const XPlane&, VertexIndex*&, VertexIndex* hint = nullptr);
        VertexIndex findCornerVertex(PosTag, VertexIndex*&);

		// manipulate
		void calcSupportingPlane();
		void calcBoundingPlane();
        void refine();

        // state
        bool isValid() const { return bIsValid; }
        void invalidate() { bIsValid = false; }
        TYPE getType() const { return TRIANGLE; }

        // octree
        int get_coords_id() const { assert(tmp_id > -1); return tmp_id; }
        void load_coords(std::deque<Real[9]>&);

    protected:
		//CGALTriangle cgalTri;
        EdgeIndex eIds[3];
        VertexIndex vIds[3];

        // normal point to inside, edge orientation incoherent with sp x bp
		XPlane bPlanes[3]; 

        int tmp_id = -1;
        bool bIsValid = true;
	};

    template <class CGALPointT>
    inline CGALPointT convertToCGALPoint(const cyPointT& pt)
    {
        return CGALPointT(pt.x, pt.y, pt.z);
    }

    template <class Kernel>
    inline typename Kernel::Triangle_3 convertToCGALTriangle(const Triangle* tri)
    {
        return typename Kernel::Triangle_3(
            convertToCGALPoint<typename Kernel::Point_3>(tri->point(0)),
            convertToCGALPoint<typename Kernel::Point_3>(tri->point(1)),
            convertToCGALPoint<typename Kernel::Point_3>(tri->point(2))
        );
    }

    // 0---------------1
    //       edge 0
    class SubPolygon : public IPolygon
    {
    public:
        SubPolygon(const Triangle* tri, uint32_t d, uint32_t i = INVALID_UINT32):
            IPolygon(i, tri->meshId()),  eIds(d), vIds(d), m_degree(d), father_(tri){}

        template <class ForwardIterator>
        void constructFromVertexList(const ForwardIterator& a, const ForwardIterator& b);
        void get_vertices_for_dumping(std::vector<VertexIndex>&) const;
        VertexIndex get_rep_vertex(EdgeIndex edgeIdx) const;
        void getAllEdges(std::vector<EdgeIndex>& output) const;
        TYPE getType() const { return SUBPOLYGON; }
        bool isValid() const { return true; }

        MyEdge& edge(int i) const;
        MyVertex& vertex(int i)const;
        EdgeIndex edgeId(int i) const { return eIds[i]; }
        VertexIndex vertexId(int i) const { return vIds[i]; }
        uint32_t degree() const { return m_degree; }

    protected:
        std::vector<EdgeIndex> eIds;
        std::vector<VertexIndex> vIds;
        ExternPtr const Triangle* father_;

        const uint32_t m_degree;
    };

    class SubPolygonWithHoles : public IPolygon
    {
    public:
        SubPolygonWithHoles(const Triangle* tri, std::vector<std::vector<VertexIndex>>& loops, uint32_t i = INVALID_UINT32 );
        ~SubPolygonWithHoles() {}
        void get_vertices_for_dumping(std::vector<VertexIndex>&) const;
        void getAllEdges(std::vector<EdgeIndex>& output) const;
        VertexIndex get_rep_vertex(EdgeIndex id) const;

        TYPE getType() const { return SUBPOLYGON_WITH_HOLES; }
        bool isValid() const { return true; }

        MyEdge& edge(int i, int j) const;
        MyVertex& vertex(int i, int j)const;
        EdgeIndex edgeId(int i, int j) const { return loops_[i].eIds[j]; }
        VertexIndex vertexId(int i, int j) const { return loops_[i].vIds[j]; }

    protected:
        // 0---------------1
        //       edge 0
        struct Loop
        {
            std::vector<EdgeIndex> eIds;
            std::vector<VertexIndex> vIds;
        };
        std::vector<Loop> loops_;
        const Triangle* father_;
    };

	class RegularMesh
	{
	public:
		typedef IPolygon FaceT;
        std::vector<std::pair<uint32_t, uint32_t>> inverseMap;

	protected:
		static  GlobalData* memmgr;
		std::vector<FaceT*>  m_faces;
		uint32_t m_id;
		bool m_bInverse = false;
        cyPointT m_center, m_scale;
        XR::BoundingBox bbox_;

	public:
		static RegularMesh* loadFromFile(const char*, uint32_t id);
		static void writeFile(RegularMesh& mesh, const char*);

        RegularMesh() :m_center(0, 0, 0), m_scale(1, 1, 1) {}
		RegularMesh(const XR::OffFile& file, uint32_t meshId); // triangle mesh
        ~RegularMesh() {}

        void clearFaces();

        // csg related
		bool& inverse() { return m_bInverse; }
		const bool& inverse() const { return m_bInverse; }
		uint32_t id() const { return m_id; }

        /// do not actually change te coordinates
        void set_transform(const cyPointT& c, const cyPointT& s) { m_center = c; m_scale = s; }

		// access
		std::vector<FaceT*>& faces() { return m_faces; }
		const std::vector<FaceT*>& faces() const  {return m_faces; }
        XR::BoundingBox bbox() const { return bbox_; }

        // geometry info
		void prepareBoolean();
	};


    /////////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION
    /////////////////////////////////////////////////////////////////////////////

    template<class ForwardIterator>
    inline void SubPolygon::constructFromVertexList(const ForwardIterator & a, const ForwardIterator & b)
    {
        ForwardIterator itr = a;
        //VertexIndex v0, v1;
        auto pMem = GlobalData::getObject();

        for (int i = 0; i < degree(); i++)
        {
            //v0 = *itr; ++itr;
            //if (itr != b)
            //    v1 = *itr;
            //else
            //    v1 = *a;

            vIds[i] = pMem->get_main_vertexId(*itr);
            ++itr;
        }

        for (int i = 0; i < degree(); i++)
        {
            eIds[i] = pMem->get_edge_id_local_scope(vIds[i], vIds[(i+1)%degree()], this);
        }
    }
}
