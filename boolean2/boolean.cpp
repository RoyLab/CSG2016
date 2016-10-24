#include "precompile.h"
#define XRWY_EXPORTS

#include <string>
#include <vector>
#include <fstream>

#include <CGAL/Bbox_3.h>
#include <CGAL/bounding_box.h>

#include <macros.h>
#include "global.h"
#include "adaptive.h"
#include "CGALext.h"

#include "boolean.h"
#include "RegularMesh.h"
#include "Octree.h"
#include "csg.h"
#include "xmemory.h"
#include "xgeometry.h"

extern "C"
{
	using namespace Boolean;

	void initContext();
	void releaseContext();

	XRWY_DLL void test(std::vector<std::string>& names, std::string& expr, const std::string& output)
	{
		const int precision = 20;
		initContext(precision);

		std::vector<RegularMesh*> meshList(names.size());
		for (int i = 0; i < names.size(); i++)
			meshList[i] = RegularMesh::loadFromFile(names[i].c_str());

		RegularMesh* result = solveCSG(expr, meshList);
		RegularMesh::writeFile(*result, output.c_str());

		SAFE_DELETE(result);
		for (auto mesh : meshList)
			SAFE_DELETE(mesh);

		releaseContext();
	}

	XRWY_DLL RegularMesh* solveCSG(const std::string& expr, std::vector<RegularMesh*>& meshes)
	{
		MemoryManager* pMem = MemoryManager::getInstance();
		XR::BoundingBox aabb(pMem->points.begin(), pMem->points.end());

		cyPointT center, scale;
		for (size_t i = 0; i < 3; i++)
		{
			center[i] = (aabb.max(i) + aabb.min(i)) / 2.0;
			scale[i] = aabb.max(i) - aabb.min(i);
		}

		for (int i = 0; i < pMem->points.size(); i++)
		{
			GS::
		}

		for (auto mesh : meshes)
		{
			mesh->transformCoords(aabb, precision);
			mesh->prepareBoolean();
		}

		CSGTree<RegularMesh>* pCsg = new CSGTree<RegularMesh>;
		pCsg->createCSGTreeFromExpr(expr, meshes.data(), meshes.size());
		pCsg->makePositiveAndLeftHeavy();

		RegularMesh* csgResult = new RegularMesh;

		Octree* pOctree = new Octree(aabb, 1e-3);
		std::vector<Octree::Node*> intersectLeaves;
		pOctree->build(meshes, &intersectLeaves);

		doIntersection(meshes, intersectLeaves);
		doClassification(pCsg, meshes, csgResult);

		csgResult->invCoords(aabb);
		return csgResult;

		SAFE_DELETE(pCsg);
		SAFE_DELETE(pOctree);
	}


	void initContext(int precision)
	{
		GS::exactinit(precision);
	}

	void releaseContext()
	{
	}

}
