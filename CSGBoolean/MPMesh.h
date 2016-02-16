#pragma once
#include "macro_util.h"
#include "csgdefs.h"
#include <OpenMesh\Core\Mesh\PolyMesh_ArrayKernelT.hh>
#include <include_fade25d\Fade_2D.h>
#include <CGAL\bounding_box.h>

namespace CSG
{
    using OpenMesh::Vec3d;
    struct ISectTriangle;
    class FeitoISectZone;

    class MyPoint : public OpenMesh::Vec3d
    {

    };

    class MyVector : public OpenMesh::Vec3d
    {

    };


    struct MyTraits : OpenMesh::DefaultTraits
    {
        //typedef OpenMesh::Vec3d Point;
        typedef MyPoint Point;
        typedef MyVector Normal;

        VertexAttributes(OpenMesh::Attributes::Color);
        FaceAttributes(OpenMesh::Attributes::Normal);
        //EdgeAttributes(OpenMesh::Attributes::Status);
    };

    class MPMesh:
        public OpenMesh::PolyMesh_ArrayKernelT<MyTraits>
    {
        typedef Bbox_3 Bbox;
        typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> Father;

        COMMON_PROPERTY(Bbox, bbox);

    public:
        MPMesh(){}
        ~MPMesh(){}

        void calcBbox()
        {
            m_bbox = CGAL::bounding_box(this->points(), vertices_end());
        }

        int  ID;
        bool bInverse;

        OpenMesh::FPropHandleT<unsigned> PointInOutTestPropHandle; // 是否在内外测试中被检测过
        OpenMesh::FPropHandleT<ISectTriangle*> SurfacePropHandle; // 是否属于相交三角形
        OpenMesh::VPropHandleT<typename Father::VertexHandle> VertexIndexPropHandle; // 在结果网格中的index
        OpenMesh::FPropHandleT<int> MarkPropHandle;
    };

    //struct Int
    //{
    //    Int():value(-1){}
    //    int value;

    //    bool operator==(int other){return value == other;}
    //    bool operator==(Int other){return value == other.value;}
    //    void operator=(int other){value = other;}
    //};

    //struct MPMesh2:public MPMeshKernel
    //{
    //    MPMesh2(const GS::BaseMesh* pMesh = nullptr);
 //       ~MPMesh2(void);

    //    MPMesh2::VertexHandle add_vertex(Vec3d& v);

    //    int  ID;
    //    bool bInverse;
    //    AABBmp BBox;

 //       const GS::BaseMesh* pOrigin;
    //    Vec3d *verticesList;
    //    OpenMesh::FPropHandleT<FeitoISectZone*> SurfacePropHandle; // 是否属于相交三角形
    //    OpenMesh::FPropHandleT<int> TopologyInfo; // 是否属于相交三角形
    //    OpenMesh::VPropHandleT<Int> VertexIndexPropHandle; // 在结果网格中的index
    //    OpenMesh::FPropHandleT<int> MarkPropHandle;
    //};

    //MPMesh* ConvertToMPMesh(const GS::BaseMesh* pMesh);
    //MPMesh2* ConvertToMPMesh2(const GS::BaseMesh* pMesh);
    //MPMesh* ConvertToMPMeshChrome(const GS::BaseMesh* pMesh);
    //inline GS::double3 Vec3dToDouble3(const Vec3d& vec){return GS::double3(vec[0], vec[1], vec[2]);}
    //inline Vec3d Double3ToVec3d(const GS::double3& vec){return Vec3d(vec[0], vec[1], vec[2]);}

    //inline void GetCorners(MPMesh* pMesh, MPMesh::FaceHandle fhandle, Vec3d*& v0, Vec3d *&v1, Vec3d *&v2)
    //{
    //    auto fvItr = pMesh->fv_begin(fhandle);
    //    v0 = &(pMesh->verticesList[fvItr->idx()]);fvItr++;
    //    v1 = &(pMesh->verticesList[fvItr->idx()]);fvItr++;
    //    v2 = &(pMesh->verticesList[fvItr->idx()]);
    //}

    //inline int TestNeighborIndex(MPMesh* pMesh, MPMesh::FaceHandle faceSeed, MPMesh::FaceHandle face)
    //{
    //    auto ffItr = pMesh->ff_begin(faceSeed);
    //    if (*ffItr == face) return 1;
    //    if (*(++ffItr) == face) return 2;
    //    if (*(++ffItr) == face) return 0;
    //    else return -1;
    //}

} // namespace CSG
