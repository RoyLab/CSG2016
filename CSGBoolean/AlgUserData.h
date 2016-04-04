#pragma once
#include <list>
#include "macroutil.h"
#include "csgdefs.h"
#include "plane_reps.h"
#include "MyMesh.h"

namespace CSG
{
    enum ContextType { CT_NONE, CT_VERTEX = 1, CT_EDGE, CT_FACET };

    template <class Refs>
    struct Context
    {
        typedef typename Refs::Face_handle FH;
        typedef typename Refs::Vertex_handle VH;
        typedef typename Refs::Halfedge_handle HEH;

        Context()
        {
            type = CT_NONE;
            meshId = -1;
            eh = nullptr;
        }

        Context(const Context& other)
        {
            type = other.type;
            meshId = other.meshId;
            switch (type)
            {
            case CT_VERTEX:
                vh = new VH(*other.vh);
                break;
            case CT_EDGE:
                eh = new HEH(*other.eh);
                break;
            case CT_FACET:
                fh = new FH(*other.fh);
                break;
            default:
                eh = nullptr;
                break;
            }
        }

        ~Context(){ release(); }

        void release()
        {
            switch (type)
            {
            case CT_VERTEX:
                SAFE_DELETE(vh);
                break;
            case CT_EDGE:
                SAFE_DELETE(eh);
                break;
            case CT_FACET:
                SAFE_DELETE(fh);
                break;
            default:
                break;
            }
        }

        ContextType type;
        int32_t     meshId;

        union
        {
            HEH* eh;
            FH* fh;
            VH* vh;
        };
    };

    struct VEntity
    {
        typedef std::vector<Context<MyMesh>> ContextList;

        VEntity()
        {
            static int count = 0;
            idx = count++;
        }

        PBPoint<K> pos;
        ContextList ctx;
        int idx = -1;

        void addContext(int meshId, MyMesh::Face_handle fh, PosTag tag)
        {
            Context<MyMesh> context;
            context.meshId = meshId;
            switch (tag)
            {
            case CSG::INNER:
                context.type = CT_FACET;
                context.fh = new MyMesh::Face_handle(fh);
                break;
            case CSG::EDGE_0:
                context.type = CT_EDGE;
                context.eh = new MyMesh::Halfedge_handle(fh->edges[0]);
                break;
            case CSG::EDGE_1:
                context.type = CT_EDGE;
                context.eh = new MyMesh::Halfedge_handle(fh->edges[1]);
                break;
            case CSG::EDGE_2:
                context.type = CT_EDGE;
                context.eh = new MyMesh::Halfedge_handle(fh->edges[2]);
                break;
            case CSG::VER_0:
                context.type = CT_VERTEX;
                context.vh = new MyMesh::Vertex_handle(fh->vertices[0]);
                break;
            case CSG::VER_1:
                context.type = CT_VERTEX;
                context.vh = new MyMesh::Vertex_handle(fh->vertices[1]);
                break;
            case CSG::VER_2:
                context.type = CT_VERTEX;
                context.vh = new MyMesh::Vertex_handle(fh->vertices[2]);
                break;
            default:
                ReportError();
                break;
            }
            ctx.push_back(context);
        }

        ContextList::iterator findInContext(int meshId)
        {
            ContextList::iterator result;
            for (result = ctx.begin(); result != ctx.end(); result++)
            {
                if (result->meshId == meshId)
                    break;
            }
            return result;
        }
        
        bool operator<(const VEntity& other)
        {
            return idx < other.idx;
        }
    };


    Relation convert(CGAL::Oriented_side side);
    Relation determineEdgeBSP(MyMesh::Halfedge_handle eh, const K::Point_3& point);
    Relation determineVertexBSP(MyMesh::Vertex_handle ctx, const K::Point_3& point);
    Relation relationOfContextNonmember(Context<MyMesh>& ctx, MyMesh::Vertex_handle vh, MyMesh::Face_handle &coins);
    Relation relationOfContext(Context<MyMesh>& ctx, MyMesh::Vertex_handle vh);
    Relation relationOfContext(Context<MyMesh>& ctx, const K::Point_3& point);

}