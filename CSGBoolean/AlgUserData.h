#pragma once
#include <list>
#include "macroutil.h"
#include "csgdefs.h"
#include "plane_reps.h"
#include "MyMesh.h"

namespace CSG
{
    struct UserEData;
    //struct LoopletTable;

    enum ContextType { CT_NONE, CT_VERTEX = 1, CT_EDGE, CT_FACET };

    struct Context
    {
        Context()
        {
            type = CT_NONE;
            meshId = -1;
            eh = nullptr;
        }

        ~Context();

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

    struct VPointer
    {
        PosTag tag = NONE;
        int idx = -1;
    };

    typedef VPointer ItstLine[2];

    typedef std::list<VEntity*>     VEntities;
    typedef VEntities::iterator     VProxy;
    typedef std::list<VProxy>       VProxies;
    typedef std::list<ItstLine>     ItstLineList;

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
        //LoopletTable*           looplets = nullptr;
        ItstLineList            isectLines;
        PBTriangle<K>*          planeRep = nullptr;
        std::vector<VProxyItr>  inVertices;

        ItstTriangle(MyMesh::Face_handle fh) :planeRep(new PBTriangle<K>(fh->)){}
    };

    struct UserFData
    {
        CGAL::Triangle_3<K>     triangle;
        TriBbox                 bbox;

        ItstTriangle*           itstTri = nullptr;
    };


    inline Context::~Context()
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
}