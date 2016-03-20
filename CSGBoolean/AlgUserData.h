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

        Context()
        {
            type = CT_NONE;
            meshId = -1;
            eh = nullptr;
        }

        ~Context()
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
            UserEData* eh;
            FH* fh;
            VH* vh;
        };
    };

    struct VEntity
    {
        std::vector<Context<MyMesh>> ctx;
    };

    struct VPointer
    {
        PosTag tag = NONE;
        int idx = -1;
    };

    typedef VPointer ItstLine[2];
    typedef std::list<ItstLine> ItstLineList;

    struct ItstTriangle
    {
        ItstLineList            isectLines;
        PBTriangle<K>*          planeRep = nullptr;
        std::vector<VProxyItr>  inVertices;

        ItstTriangle(MyMesh::Face_handle fh){}
    };

}