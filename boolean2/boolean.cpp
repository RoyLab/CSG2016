#include "precompile.h"
#include <string>
#include <vector>
#include <fstream>

#define XRWY_EXPORTS
#include <macros.h>
#include "global.h"
#include "adaptive.h"
#include "RegularSetMesh.h"
#include "csg.h"

extern "C"
{
	using namespace Boolean;
	using namespace CSG;

	void initContext();
	void releaseContext();

	void XRWY_DLL test(std::vector<std::string>& names, std::string& expr, const std::string& output)
	{
		initContext();

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

	RegularMesh* XRWY_DLL solveCSG(const std::string& expr, std::vector<RegularMesh*>& meshes)
	{
		pMeshList = &meshes;
		for (auto mesh : *pMeshList)
			mesh->calcBbox();

		CGAL::Bbox_3 aabb = pMeshList->at(0)->get_bbox();
		for (auto mesh : *pMeshList)
			aabb += mesh->get_bbox();

		for (auto mesh : *pMeshList)
		{
			mesh->normalize(aabb);
			mesh->filter();
			mesh->init();
		}

		CSGTree<MyMesh>* pCsg = new CSGTree<MyMesh>;
		pCsg->createCSGTreeFromExpr(expr, pMeshList->data(), pMeshList->size());
		pCsg->makePositiveAndLeftHeavy();

		Octree* pOctree = new Octree;
		std::vector<Octree::Node*> intersectLeaves;
		pOctree->build(*pMeshList, &intersectLeaves);

		itst = new ItstAlg(pMeshList);
		itst->doIntersection(intersectLeaves);

		//std::cout.precision(18);
		//for (auto pt : itst->vEnt)
		//std::cout << "point: " << pt->pos.getCoord() << std::endl;

		floodColoring(pCsg, itst);

		csgResult->denormalize(aabb);

		SAFE_DELETE(itst);
		SAFE_DELETE(pCsg);
		SAFE_DELETE(pOctree);
	}


	void initContext()
	{
		GS::exactinit();
	}

	void releaseContext()
	{
	}

}
