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

        bool hasContext(size_t meshId) const
        {
            return findInContext(meshId) != ctx.end();
        }

        ContextList::const_iterator findInContext(size_t meshId) const 
        {
            ContextList::const_iterator result;
            for (result = ctx.cbegin(); result != ctx.cend(); result++)
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

    Relation relationOfContext(const Context<MyMesh>& ctx, const PBPoint<K>& point, Context<MyMesh>** pCtx = nullptr);
    Relation determineRelationOfFacet(const Context<MyMesh>& ctx, const PBPoint<K>& p0, const PBPoint<K>& p1);
}