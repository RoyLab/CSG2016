#include "AlgUserData.h"
#include "BSPTree.h"
#include "MyMesh.h"

namespace CSG
{
    DECLARE_MYMESH_TYPES;


    Relation convert(CGAL::Oriented_side side)
    {
        switch (side)
        {
        case CGAL::ON_NEGATIVE_SIDE:
            return REL_INSIDE;
        case CGAL::ON_POSITIVE_SIDE:
            return REL_OUTSIDE;
        case CGAL::ON_ORIENTED_BOUNDARY:
            return REL_ON_BOUNDARY;
        default:
            ReportError();
            return REL_NOT_AVAILABLE;
        }
    }

    /* canb be coplanar */
    Relation determineFaceBSP(FH fh, const PBPoint<K>& point, Context<MyMesh>** pCtx)
    {
        auto side = fh->data->sp.oriented_side(point.getCoord());
        auto result = convert(side);

        if (result == REL_ON_BOUNDARY && pCtx)
        {
            (*pCtx)->type = CT_FACET;
            (*pCtx)->fh = new FH(fh);
        }
        return result;
    }


    Relation determineEdgeBSP(EH eh, const PBPoint<K>& point, Context<MyMesh>** pCtx)
    {
        auto f0 = eh->facet();
        auto f1 = eh->opposite()->facet();

        CGAL::Sign r0 = f0->data->sp.oriented_side(point.getCoord());
        CGAL::Sign r1 = f1->data->sp.oriented_side(point.getCoord());

        Relation result;

        if (r0 == r1) result = convert(r0);

        if (r0 == CGAL::ON_ORIENTED_BOUNDARY)
        {
            auto opvh = eh->next()->vertex();
            CGAL::Sign r01 = f1->data->sp.oriented_side(opvh->point());
            if (r01 == r1)
                result = convert(r0);
            else
                result = convert(r1);
        }
        else
        {
            auto opvh = eh->opposite()->next()->vertex();
            CGAL::Sign r10 = f0->data->sp.oriented_side(opvh->point());
            if (r10 == r0)
                result = convert(r1);
            else
                result = convert(r0);
        }

        if (pCtx && result == REL_ON_BOUNDARY)
        {
            auto ptr = *pCtx;
            assert(ptr);
            if (r0 == r1)
            {
                ptr->type = CT_EDGE;
                ptr->eh = new EH(eh);
            }
            else
            {
                ptr->type = CT_FACET;
                if (r0 == CGAL::ON_ORIENTED_BOUNDARY)
                    ptr->fh = new FH(f0);
                else
                    ptr->fh = new FH(f1);
            }
        }

        return result;
    }

    /* not implemented yet! */
    Relation determineVertexBSP(VH ctx, const PBPoint<K>& point, Context<MyMesh>** pCtx)
    {
        ReportError("Not implemented vertex context classfication!");
        std::vector<FH> facets;
        auto end = ctx->halfedge();
        auto cur = end;
        do
        {
            facets.push_back(cur->facet());
            cur = cur->next()->opposite();
        } while (cur != end);

        BSPTree* pTree = new BSPTree(facets);
        SAFE_DELETE(pTree);

        return REL_NOT_AVAILABLE;
    }

    Relation relationOfContext(const Context<MyMesh>& ctx, const PBPoint<K>& point, Context<MyMesh>** pCtx)
    {
        Relation result = REL_NOT_AVAILABLE;
        if (pCtx && *pCtx)
            (*pCtx)->meshId = ctx.meshId;

        switch (ctx.type)
        {
        case CT_VERTEX:
            result = determineVertexBSP(*ctx.vh, point, pCtx);
            break;
        case CT_EDGE:
            result = determineEdgeBSP(*ctx.eh, point, pCtx);
            break;
        case CT_FACET:
            result = determineFaceBSP(*ctx.fh, point, pCtx);
            break;
        default:
            ReportError();
            break;
        }
        return result;
    }


    Relation determineRelationOfFacet(const Context<MyMesh>& ctx, const PBPoint<K>& p0, const PBPoint<K>& p1)
    {
        Context<MyMesh> *eCtx = new Context<MyMesh>;
        Relation rel = relationOfContext(ctx, p0, &eCtx);

        if (rel == REL_ON_BOUNDARY)
            rel = relationOfContext(*eCtx, p1);

        SAFE_DELETE(eCtx);
        return rel;
    }


}