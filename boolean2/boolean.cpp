#include "precompile.h"
#include <string>
#include <vector>
#include <fstream>

#include <CGAL/Bbox_3.h>

#define XRWY_EXPORTS
#include <macros.h>
#include "global.h"
#include "adaptive.h"
#include "CGALext.h"

#include "boolean.h"
#include "RegularSetMesh.h"
#include "csg.h"

extern "C"
{
	using namespace Boolean;

	void initContext();
	void releaseContext();

	XRWY_DLL void test(std::vector<std::string>& names, std::string& expr, const std::string& output)
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

	XRWY_DLL RegularMesh* solveCSG(const std::string& expr, std::vector<RegularMesh*>& meshes)
	{
		CGAL::Bbox_3 aabb = meshes[0]->bbox();
		for (auto meshItr = meshes.begin()+1; 
			meshItr != meshes.end(); meshItr++)
			aabb += (*meshItr)->bbox();

		for (auto mesh : meshes)
		{
			mesh->transformCoords(aabb, 20);
			mesh->init();
		}

		CSGTree<RegularMesh>* pCsg = new CSGTree<RegularMesh>;
		pCsg->createCSGTreeFromExpr(expr, meshes.data(), meshes.size());
		pCsg->makePositiveAndLeftHeavy();

		RegularMesh* csgResult = new RegularMesh;

		Octree* pOctree = new Octree;
		std::vector<Octree::Node*> intersectLeaves;
		pOctree->build(*pMeshList, &intersectLeaves);

		//itst = new ItstAlg(pMeshList);
		//itst->doIntersection(intersectLeaves);

		////std::cout.precision(18);
		////for (auto pt : itst->vEnt)
		////std::cout << "point: " << pt->pos.getCoord() << std::endl;

		//floodColoring(pCsg, itst);

		csgResult->invCoords(aabb);
		return csgResult;

		//SAFE_DELETE(itst);
		//SAFE_DELETE(pCsg);
		//SAFE_DELETE(pOctree);
	}


	void initContext()
	{
		GS::exactinit();
	}

	void releaseContext()
	{
	}

}
