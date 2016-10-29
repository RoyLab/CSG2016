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
#include "intersection.h"

extern "C"
{
	using namespace Boolean;

	void initContext(int precision);
	void releaseContext();

	XRWY_DLL void test(std::vector<std::string>& names, std::string& expr, const std::string& output)
	{
		const int precision = 16;
		initContext(precision);

        test1();
        test2();
        test3();
        test4();

		std::vector<RegularMesh*> meshList(names.size());
		for (int i = 0; i < names.size(); i++)
		{
			meshList[i] = RegularMesh::loadFromFile(names[i].c_str(), i);
		}

		RegularMesh* result = solveCSG(expr, meshList);
		//RegularMesh::writeFile(*result, output.c_str());

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
		for (int i = 0; i < 3; i++)
		{
			center[i] = (aabb.max(i) + aabb.min(i)) / 2.0;
			scale[i] = aabb.max(i) - aabb.min(i);
		}

		for (auto &pt: pMem->points)
		{
			XR::normalizeCoords(center, scale, pt);
			fp_filter(reinterpret_cast<Real*>(&pt));
		}

		for (auto mesh : meshes)
			mesh->prepareBoolean();

		CSGTree<RegularMesh>* pCsg = new CSGTree<RegularMesh>;
		pCsg->createCSGTreeFromExpr(expr, meshes.data(), meshes.size());
		pCsg->makePositiveAndLeftHeavy();

		RegularMesh* csgResult = new RegularMesh;

		Octree* pOctree = new Octree();
		const double padding = 1e-3;
		std::vector<Octree::Node*> intersectLeaves;

		auto cgalbbox = Bbox_3(-1, -1, -1, 1, 1, 1);
		cgalbbox = enlarge(cgalbbox, padding);
		pOctree->build(meshes, cgalbbox, &intersectLeaves);

		doIntersection(meshes, intersectLeaves);
		//doClassification(pCsg, meshes, csgResult);

		//csgResult->invCoords(aabb);
		//return csgResult;

		//SAFE_DELETE(pCsg);
		//SAFE_DELETE(pOctree);

		return csgResult;
	}


	void initContext(int precision)
	{
		exactinit(precision);
	}

	void releaseContext()
	{
	}

}
