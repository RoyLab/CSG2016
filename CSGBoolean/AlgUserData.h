#pragma once
#include <list>
#include "macroutil.h"
#include "csgdefs.h"
#include "MyMesh.h"

namespace CSG
{
    struct UserEData;
    struct LoopletTable;

    enum ContextType { CT_NONE, CT_VERTEX = 1, CT_EDGE, CT_FACET };

    struct Context
    {
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
            MyMesh::Face_handle* fh;
            MyMesh::Vertex_handle* vh;
        };
    };

    struct VEntity
    {
        std::vector<Context> ctx;
    };

    typedef std::list<VEntity*>     VEntities;
    typedef VEntities::iterator     VProxy;
    typedef std::list<VProxy>       VProxies;

    class VProxyItr :
        public std::list<VProxy>::iterator
    {
    public:

        VEntity* pointer()
        {
            return ***this;
        }

        const VEntity* pointer() const
        {
            return ***this;
        }
    };

    typedef Cube_3 TriBbox;

    struct UserVData
    {
        VProxyItr* proxy = nullptr;
    };

    struct UserEData
    {
        std::vector<VProxyItr>  vertices;
    };

    struct ItstTriangle
    {
        LoopletTable*           looplets = nullptr;
        std::vector<VProxyItr>  inVertices;
    };

    struct UserFData
    {
        CGAL::Triangle_3<K>     triangle;
        TriBbox                 bbox;

        ItstTriangle*           itstTri = nullptr;
    };
}