#include "precompile.h"
#define XRWY_EXPORTS

#include <string>
#include <vector>
#include <fstream>

#include <CGAL/Bbox_3.h>
#include <CGAL/bounding_box.h>

#include <macros.h>
#include <XTimer.hpp>
#include <xlogger.h>
#include "global.h"
#include "adaptive.h"
#include "CGALext.h"

#include "boolean.h"
#include "RegularMesh.h"
#include "Octree.h"
#include "csg.h"
#include "intersection.h"


namespace Boolean
{
    void tessellation(std::vector<Boolean::RegularMesh*>& meshes);
    void doClassification(Octree*, CSGTree<RegularMesh>*, std::vector<RegularMesh*>&, RegularMesh*, MyVertex::Index);

    void initContext() {}
    void releaseContext() {}

    bool collinear(const cyPointT& a, const cyPointT& b, const cyPointT& c)
    {
        cyPointT e[2] = { c - a, c - b };
        cyPointT n = e[0].Cross(e[1]);
        return n.x == 0 && n.y == 0 && n.z == 0;
    }

    MyVertex::Index pickSeed(std::vector<MyVertex::Index>& seed)
    {
        cyPointT basePoint[2];
        for (uint32_t i = 0; i < seed.size(); i++)
        {
            MyVertex& vRef = xvertex(seed[i]);
            basePoint[0] = vRef.point();
            auto eItr = vRef.edges.begin();

            MyVertex::Index base[2];
            MyEdge& edge = xedge(*eItr++);
            MyVertex::Index theother =
                edge.ends[0] == seed[i] ? edge.ends[1] : edge.ends[0];
            base[0] = theother;
            basePoint[1] = xvertex(base[0]).point();

            while (1)
            {
                MyEdge& edge = xedge(*eItr++);
                MyVertex::Index theother =
                    edge.ends[0] == seed[i] ? edge.ends[1] : edge.ends[0];

                if (theother != base[0] && !collinear(basePoint[0], basePoint[1], xvertex(theother).point()))
                {
                    base[1] = theother;
                    break;
                }
            }

            assert(!xvertex(base[0]).isPlaneRep());
            assert(!xvertex(base[1]).isPlaneRep());
            assert(eItr != vRef.edges.end());

            XPlaneBase pbase(basePoint[0], basePoint[1], xvertex(base[1]).point());
            bool got = false;
            while (eItr != vRef.edges.end())
            {
                MyEdge& edge = xedge(*eItr++);
                MyVertex::Index theother =
                    edge.ends[0] == seed[i] ? edge.ends[1] : edge.ends[0];
                MyVertex& vRef2 = xvertex(theother);
                if (!vRef2.isPlaneRep() && pbase.orientation(vRef2.point()) != ON_ORIENTED_BOUNDARY)
                {
                    got = true;
                    break;
                }
            }

            if (got)
                return seed[i];
        }
        assert(0);
        return INVALID_UINT32;
    }
}


extern "C"
{
	using namespace Boolean;

	XRWY_DLL void test(std::vector<std::string>& names, std::string& expr, const std::string& output)
	{
        const double dieta = 1;
        const int precision = exactinit(dieta);
        initContext();

        test1();
        test2();
        test3();
        test4();

		std::vector<RegularMesh*> meshList(names.size());
		for (int i = 0; i < names.size(); i++)
		{
			meshList[i] = RegularMesh::loadFromFile(names[i].c_str(), i);
		}

        XTIMER_HELPER(setClock("main"));
		RegularMesh* result = solveCSG(expr, meshList);
        XLOG_INFO << "Overall time: " << XTIMER_HELPER(milliseconds("main")) << " ms";

		RegularMesh::writeFile(*result, output.c_str());

		SAFE_DELETE(result);
		for (auto mesh : meshList)
			SAFE_DELETE(mesh);

		releaseContext();
	}

	XRWY_DLL RegularMesh* solveCSG(const std::string& expr, std::vector<RegularMesh*>& meshes)
	{
		MemoryManager* pMem = MemoryManager::getInstance();
        std::vector<MyVertex::Index> xmins;
		XR::BoundingBox aabb(pMem->points.begin(), pMem->points.end(), xmins);

		cyPointT center, scale;
		for (int i = 0; i < 3; i++)
		{
			center[i] = (aabb.maxVal(i) + aabb.minVal(i)) / 2.0;
			scale[i] = aabb.maxVal(i) - aabb.minVal(i);
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

		auto cgalbbox = Bbox_3(-1,-1,-1, 1, 1, 1);
		cgalbbox = enlarge(cgalbbox, padding);
		pOctree->build(meshes, cgalbbox, &intersectLeaves);

		doIntersection(meshes, intersectLeaves);
        MyVertex::Index seed = pickSeed(xmins);
        //pMem->outputIntersection("C:/Users/XRwy/Desktop/x2.xyz", center, scale);

        tessellation(meshes);

        //meshes[0]->invCoords(center, scale);
        //RegularMesh::writeFile(*meshes[0], "D:/a.off");

        doClassification(pOctree, pCsg, meshes, csgResult, seed);
		csgResult->invCoords(center, scale);

		SAFE_DELETE(pCsg);
		SAFE_DELETE(pOctree);

		return csgResult;
	}



}
