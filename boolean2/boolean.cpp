#include "precompile.h"
#define XRWY_EXPORTS
#define XTIMER_INSTANCE
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

//#include <vld.h>

namespace Boolean
{
    extern int xrcount, pre_count, after_count;
    void tessellation(std::vector<Boolean::RegularMesh*>& meshes, std::vector<Triangle*> insct_triangles);
    void doClassification(Octree*, CSGTree<RegularMesh>*, std::vector<RegularMesh*>&, RegularMesh*, VertexIndex, RegularMesh*);
    void doIntersection(std::vector<RegularMesh*>&, std::vector<Octree::Node*>&, std::vector<Triangle*>& insct_triangles);

    void initContext() {}
    void releaseContext()
    {
        GlobalData::getObject()->clear();
    }

    bool collinear(const cyPointT& a, const cyPointT& b, const cyPointT& c)
    {
        cyPointT e[2] = { c - a, c - b };
        cyPointT n = e[0].Cross(e[1]);
        return n.x == 0 && n.y == 0 && n.z == 0;
    }

    // 取一个点，它的一度点不都在一个平面上
    VertexIndex pickSeed(std::vector<VertexIndex>& seed)
    {
        cyPointT basePoint[2];
        for (uint32_t i = 0; i < seed.size(); i++)
        {
            if (GlobalData::getObject()->get_main_vertexId(seed[i]) != seed[i])
            {
                continue;
            }

            MyVertex& vRef = xvertex(seed[i]);
            basePoint[0] = vRef.vertex_rep();
            auto eItr = vRef.edges_local().begin();

            VertexIndex base[2];
            MyEdge& edge = xedge(*eItr++);
            VertexIndex theother = edge.theOtherVId(seed[i]);
                //edge.ends[0] == seed[i] ? edge.ends[1] : edge.ends[0];
            base[0] = theother;
            basePoint[1] = xvertex(base[0]).vertex_rep();

            while (1)
            {
                MyEdge& edge = xedge(*eItr++);
                VertexIndex theother = edge.theOtherVId(seed[i]);
                    //edge.ends[0] == seed[i] ? edge.ends[1] : edge.ends[0];

                if (!vertex_id_equals_simple(theother, base[0]) && !collinear(basePoint[0], basePoint[1], xvertex(theother).vertex_rep()))
                {
                    base[1] = theother;
                    break;
                }
            }

            assert(!xvertex(base[0]).isPlaneRep());
            assert(!xvertex(base[1]).isPlaneRep());
            assert(eItr != vRef.edges_local().end());

            XPlane pbase;
            pbase.construct_from_three_vertices(basePoint[0], basePoint[1], xvertex(base[1]).vertex_rep());
            bool got = false;
            while (eItr != vRef.edges_local().end())
            {
                MyEdge& edge = xedge(*eItr++);
                VertexIndex theother = edge.theOtherVId(seed[i]);
                    //edge.ends[0] == seed[i] ? edge.ends[1] : edge.ends[0];
                MyVertex& vRef2 = xvertex(theother);
                if (!vRef2.isPlaneRep() && pbase.orientation(vRef2.vertex_rep()) != ON_ORIENTED_BOUNDARY)
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

	XRWY_DLL void test(const CsgInputs& inputs)
	{
        XLOG_INFO << "\n***Test Begin***";
        XLOG_INFO << inputs.names.size() << " meshes:";
        for (int i = 0; i < inputs.names.size(); i++)
            XLOG_INFO << inputs.names[i];

#ifdef XR_PROFILE
        XLOG_INFO << "______tHIS iS pROFILING vERTION_____";
#endif
        GlobalData* pMem = GlobalData::getObject();
        initContext();

        auto& meshlist = pMem->meshes;
        meshlist.resize(inputs.names.size());
		for (int i = 0; i < inputs.names.size(); i++)
		{
            meshlist[i] = RegularMesh::loadFromFile(inputs.names[i].c_str(), i);
		}

        if (inputs.need_pwn_test)
        {
            pMem->results.is_pwn = meshlist[0]->is_pwn();
        }

        XTIMER_HELPER(setClock("main"));
        RegularMesh* result = nullptr;
        try
        {
            result = solveCSG(inputs.expr, meshlist);
            pMem->results.total_time = XTIMER_HELPER(milliseconds("main"));
            XLOG_INFO << "Overall time: " << pMem->results.total_time << " ms";
            RegularMesh::writeFile(*result, inputs.output.c_str());
            pMem->results.error = 0;
        }
        catch (...)
        {
            XLOG_FATAL << "Unexpected error.";
            pMem->results.n_faces = 0;
            pMem->results.n_vertices = 0;
            pMem->results.error = 1;
        }

		SAFE_DELETE(result);
        for (auto mesh : meshlist)
        {
            mesh->clearFaces();
			SAFE_DELETE(mesh);
        }

        // record results
        if (inputs.need_log)
        {
            std::ofstream log_file(inputs.log);
            log_file << pMem->results.total_time << std::endl;
            log_file << pMem->results.n_faces << std::endl;
            log_file << pMem->results.n_vertices << std::endl;
            log_file << pMem->results.error << std::endl;
            for (int i = 0; i < 4; ++i)
            {
                log_file << pMem->results.step_time[i] << std::endl;
            }

            if (inputs.need_pwn_test)
            {
                log_file << pMem->results.is_pwn << std::endl;
            }

            log_file.close();
        }

		releaseContext();
	}

    XRWY_DLL bool is_pwn(const char * name)
    {
        RegularMesh *mesh = RegularMesh::loadFromFile(name, 0);
        if (!mesh)
        {
            XLOG_FATAL << "Cannot load the mesh " << name;
            throw std::exception();
        }

        bool ret = mesh->is_pwn();
        delete mesh;
        return ret;
    }

	XRWY_DLL RegularMesh* solveCSG(const std::string& expr, std::vector<RegularMesh*>& meshes)
	{
        if (meshes.size() > MAX_MESH_COUNT)
        {
            XLOG_FATAL << "Mesh number exceeds limit";
            exit(-1);
        }

		GlobalData* pMem = GlobalData::getObject();
        XR::BoundingBox aabb(meshes[0]->bbox());
        for (int i = 1; i < meshes.size(); ++i)
        {
            aabb.include_bbox(meshes[i]->bbox());
        }

        std::vector<VertexIndex> xmins;
        double xmin = 1e99;
        for (int i = 0; i < pMem->vertices.size(); ++i)
        {
            double x = pMem->vertices[i].vertex_rep().x;
            if (x <= xmin)
            {
                if (xmin > x)
                {
                    xmin = x;
                    xmins.clear();
                }
                xmins.push_back(i);
            }
        }

        double logDieta = -2;
        if (pMem->points.size() > 1e6) logDieta = -5;
        else if (pMem->points.size() < 1000) logDieta = 0;
        exactinit(logDieta);

        //std::vector<VertexIndex> xmins2;
        //XR::BoundingBox aabb2(pMem->points.begin(), pMem->points.end(), xmins2);

        // transform vertices coordinates
		cyPointT center, scale;
		for (int i = 0; i < 3; i++)
		{
			center[i] = (aabb.maxVal(i) + aabb.minVal(i)) / 2.0;
			scale[i] = aabb.maxVal(i) - aabb.minVal(i);
		}

		for (int i = 0; i < pMem->points.size(); i++)
		{
			XR::normalizeCoords(center, scale, pMem->points[i]);
			fp_filter(reinterpret_cast<Real*>(&pMem->points[i]));
		}

        for (int i = 0; i < meshes.size(); ++i)
        {
            meshes[i]->set_transform(center, scale);
        }

        // do some initialization
		for (auto mesh : meshes)
			mesh->prepareBoolean();

		CSGTree<RegularMesh>* pCsg = new CSGTree<RegularMesh>;
		pCsg->createCSGTreeFromExpr(expr, meshes.data(), meshes.size());
		pCsg->makePositiveAndLeftHeavy();

        XTIMER_HELPER(setClock("octree"));
        Octree* pOctree = new Octree();
		const double padding = 1e-3;
		std::vector<Octree::Node*> intersectLeaves;

		auto cgalbbox = Bbox_3(-1,-1,-1, 1, 1, 1);
		cgalbbox = enlarge(cgalbbox, padding);
		pOctree->build(meshes, cgalbbox, false, &intersectLeaves);
        pMem->results.step_time[0] = XTIMER_HELPER(milliseconds("octree"));
        XLOG_INFO << "build octree time: " << pMem->results.step_time[0] << " ms";
        XLOG_INFO << "Number of triaabb test: " << xrcount <<" / " << pre_count << " / " << after_count;


        XTIMER_HELPER(setClock("insct"));
        std::vector<Triangle*> insct_triangles;
        doIntersection(meshes, intersectLeaves, insct_triangles);
        pMem->results.step_time[1] = XTIMER_HELPER(milliseconds("insct"));
        XLOG_INFO << "intersection test time: " << pMem->results.step_time[1] << " ms";

        XTIMER_HELPER(setClock("tess"));
        tessellation(meshes, insct_triangles);
        pMem->results.step_time[2] = XTIMER_HELPER(milliseconds("tess"));
        XLOG_INFO << "tessellation time: " << pMem->results.step_time[2] << " ms";

        //pMem->outputIntersection("C:/Users/XRwy/Desktop/x2.xyz", center, scale);
        //meshes[0]->invCoords(center, scale);

        //for (int i = 0; i < 1; ++i)
        //{
        //    char a[80];
        //    sprintf(a, "F:/a%d.off", i);

        //    RegularMesh::writeFile(*meshes[i], a);
        //}
        //exit(0);

        XTIMER_HELPER(setClock("classify"));
        VertexIndex seed = pickSeed(xmins);
        RegularMesh* csgResult = new RegularMesh;
        RegularMesh* debug_mesh = new RegularMesh;

        doClassification(pOctree, pCsg, meshes, csgResult, seed, debug_mesh);

        pMem->results.step_time[3] = XTIMER_HELPER(milliseconds("classify"));
        XLOG_INFO << "classification time: " << pMem->results.step_time[3] << " ms";

        csgResult->set_transform(center, scale);
        debug_mesh->set_transform(center, scale);

        //RegularMesh::writeFile(*debug_mesh, "D:/debug.off");
        delete debug_mesh;

		SAFE_DELETE(pCsg);
		SAFE_DELETE(pOctree);

		return csgResult;
	}



}
