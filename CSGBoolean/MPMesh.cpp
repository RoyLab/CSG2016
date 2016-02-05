#include "precompile.h"
#include "MPMesh.h"
#include <limits>
#include "BaseMesh.h"
#include "COctree.h"
#include "Bool.h"

namespace CSG
{
	static const unsigned VERTEX_ID_INIT = UINT_MAX;

    extern "C" CSG_API void ReleaseMPMesh2(MPMesh2*& ptr) {SAFE_RELEASE(ptr);}

	inline OpenMesh::Vec3d convert_double3(GS::double3& vec)
	{
		return OpenMesh::Vec3d(vec.x, vec.y, vec.z);
	}

	MPMesh2::MPMesh2(const GS::BaseMesh* pMesh):
        ID(-1), pOrigin(pMesh), bInverse(false),
		verticesList(nullptr)
    {
		request_face_normals();
        add_property(MarkPropHandle);
        add_property(SurfacePropHandle);
		add_property(VertexIndexPropHandle);
		add_property(TopologyInfo);
    }

    MPMesh2::~MPMesh2(void)
    {
		release_face_normals();
        remove_property(MarkPropHandle);
        remove_property(SurfacePropHandle);
		remove_property(VertexIndexPropHandle);
		remove_property(TopologyInfo);
		SAFE_RELEASE_ARRAY(verticesList);
    }


	MPMesh2::VertexHandle MPMesh2::add_vertex(Vec3d& v)
	{
		BBox.IncludePoint(v);
		return MPMeshKernel::add_vertex(v);
	}

	MPMesh::MPMesh(const GS::BaseMesh* pMesh):
        ID(-1), pOrigin(pMesh), bInverse(false),
		verticesList(nullptr)
    {
		request_face_normals();
		add_property(PointInOutTestPropHandle);
		add_property(SurfacePropHandle);
		add_property(VertexIndexPropHandle);
		add_property(MarkPropHandle);
    }

    MPMesh::~MPMesh(void)
    {
		release_face_normals();
		remove_property(PointInOutTestPropHandle);
		remove_property(SurfacePropHandle);
		remove_property(VertexIndexPropHandle);
		remove_property(MarkPropHandle);
		SAFE_RELEASE_ARRAY(verticesList);
    }

	MPMesh* ConvertToMPMeshChrome(const GS::BaseMesh* pMesh)
	{
		MPMesh *res = new MPMesh(pMesh);
		res->request_face_colors();
		int n = (int)pMesh->VertexCount();
		res->BBox.Clear();
		MPMesh::VertexHandle vhandle;
		res->verticesList = new Vec3d[n];
		for (int i = 0; i < n; i++)
		{
			res->verticesList[i] = convert_double3(pMesh->Vertex(i).pos);
			vhandle = res->add_vertex(res->verticesList[i]);
			//res->property(res->PointInOutTestPropHandle, vhandle) = 0;
			res->BBox.IncludePoint(res->point(vhandle));
		}

		std::vector<MPMesh::VertexHandle>  face_vhandles(3);
		n = (int)pMesh->PrimitiveCount();
		for (int i = 0; i < n; i++)
		{
			auto &info = pMesh->TriangleInfo(i);
			face_vhandles[0] = res->vertex_handle(info.VertexId[0]);
			face_vhandles[1] = res->vertex_handle(info.VertexId[1]);
			face_vhandles[2] = res->vertex_handle(info.VertexId[2]);
			MPMesh::FaceHandle fh = res->add_face(face_vhandles);
			res->set_color(fh, MPMesh::Color(info.color[0], info.color[1], info.color[2]));
		}

		res->update_normals();
		return res;
	}


	MPMesh* ConvertToMPMesh(const GS::BaseMesh* pMesh)
	{
		MPMesh *res = new MPMesh(pMesh);

        int n = (int)pMesh->VertexCount();
        res->BBox.Clear();
		MPMesh::VertexHandle vhandle;
		res->verticesList = new Vec3d[n];
        for (int i = 0; i < n; i++)
		{
			res->verticesList[i] = convert_double3(pMesh->Vertex(i).pos);
            vhandle = res->add_vertex(res->verticesList[i]);
			//res->property(res->PointInOutTestPropHandle, vhandle) = 0;
			res->BBox.IncludePoint(res->point(vhandle));
		}

		std::vector<MPMesh::VertexHandle>  face_vhandles(3);
        n = (int)pMesh->PrimitiveCount();
        for (int i = 0; i < n; i++)
        {
            auto &info = pMesh->TriangleInfo(i);
			face_vhandles[0] = res->vertex_handle(info.VertexId[0]);
			face_vhandles[1] = res->vertex_handle(info.VertexId[1]);
			face_vhandles[2] = res->vertex_handle(info.VertexId[2]);
            res->add_face(face_vhandles);
        }

		res->update_normals();

		return res;
	}


    	MPMesh2* ConvertToMPMesh2(const GS::BaseMesh* pMesh)
	{
		MPMesh2 *res = new MPMesh2(pMesh);

        int n = (int)pMesh->VertexCount();
        res->BBox.Clear();
		MPMesh2::VertexHandle vhandle;
		res->verticesList = new Vec3d[n];
        for (int i = 0; i < n; i++)
		{
			res->verticesList[i] = convert_double3(pMesh->Vertex(i).pos);
            vhandle = res->add_vertex(res->verticesList[i]);
			res->BBox.IncludePoint(res->point(vhandle));
		}

		std::vector<MPMesh2::VertexHandle>  face_vhandles(3);
        n = (int)pMesh->PrimitiveCount();
        for (int i = 0; i < n; i++)
        {
            auto &info = pMesh->TriangleInfo(i);
			face_vhandles[0] = res->vertex_handle(info.VertexId[0]);
			face_vhandles[1] = res->vertex_handle(info.VertexId[1]);
			face_vhandles[2] = res->vertex_handle(info.VertexId[2]);
            res->add_face(face_vhandles);
        }

		res->update_normals();

		return res;
	}
} // namespace CSG 
