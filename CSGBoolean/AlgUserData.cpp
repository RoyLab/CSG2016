#include "AlgUserData.h"
#include "BSPTree.h"

namespace CSG
{
    Relation convert(CGAL::Oriented_side side)
    {
        switch (side)
        {
        case CGAL::ON_NEGATIVE_SIDE:
            return REL_INSIDE;
        case CGAL::ON_POSITIVE_SIDE:
            return REL_OUTSIDE;
        case CGAL::ON_ORIENTED_BOUNDARY:
            return REL_SAME;
        default:
            ReportError();
            return REL_NOT_AVAILABLE;
        }
    }

    Relation determineEdgeBSP(MyMesh::Halfedge_handle eh, const K::Point_3& point)
    {
        auto f0 = eh->facet();
        auto f1 = eh->opposite()->facet();

        auto r0 = f0->data->sp.oriented_side(point);
        auto r1 = f1->data->sp.oriented_side(point);

        if (r0 == r1) return convert(r0);

        if (r0 == REL_SAME)
        {
            auto opvh = eh->next()->vertex();
            auto r01 = f1->data->sp.oriented_side(opvh->point());
            if (r01 == r1)
                return convert(r0);
            else
                return convert(r1);
        }
        else
        {
            auto opvh = eh->opposite()->next()->vertex();
            auto r10 = f0->data->sp.oriented_side(opvh->point());
            if (r10 == r0)
                return convert(r1);
            else
                return convert(r0);
        }
    }

    Relation determineVertexBSP(MyMesh::Vertex_handle ctx, const K::Point_3& point)
    {
        typedef MyMesh::Face_handle FH;
        std::vector<FH> facets;
        auto end = ctx->halfedge();
        auto cur = end;
        do
        {
            facets.push_back(cur->facet());
            cur = cur->next()->opposite();
        } while (cur != end);

        BSPTree* pTree = new BSPTree(facets);
        ReportError("Not implemented!");
        SAFE_DELETE(pTree);

        return REL_NOT_AVAILABLE;
    }

    Relation relationOfContextNonmember(Context<MyMesh>& ctx, MyMesh::Vertex_handle vh, MyMesh::Face_handle &coins)
    {
        CGAL::Oriented_side side, side1;
        Relation result = REL_NOT_AVAILABLE;
        switch (ctx.type)
        {
        case CT_VERTEX:
            result = determineVertexBSP(*ctx.vh, vh->point());
        case CT_EDGE:
            result = determineEdgeBSP(*ctx.eh, vh->point());
            break;
        case CT_FACET:
            side = (*ctx.fh)->data->sp.oriented_side(vh->point());
            result = convert(side);
            if (result == REL_SAME) coins = (*ctx.fh);
            break;
        default:
            ReportError();
            break;
        }
        return result;
    }

    Relation relationOfContext(Context<MyMesh>& ctx, MyMesh::Vertex_handle vh)
    {
        MyMesh::Face_handle coins;
        return relationOfContextNonmember(ctx, vh, coins);
    }

    Relation relationOfContext(Context<MyMesh>& ctx, const K::Point_3& point)
    {
        CGAL::Oriented_side side, side1;
        Relation result = REL_NOT_AVAILABLE;
        switch (ctx.type)
        {
        case CT_VERTEX:
            result = determineVertexBSP(*ctx.vh, point);
        case CT_EDGE:
            result = determineEdgeBSP(*ctx.eh, point);
            break;
        case CT_FACET:
            side = (*ctx.fh)->data->sp.oriented_side(point);
            result = convert(side);
            break;
        default:
            ReportError();
            break;
        }
        return result;
    }

}