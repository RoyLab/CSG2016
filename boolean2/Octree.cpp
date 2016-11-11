#include "precompile.h"
#include "Octree.h"
#include "CGALext.h"
#include "RegularMesh.h"
#include "xgeometry.h"
#include "triboxtest.hpp"


namespace Boolean
{
    int xrcount;
    namespace
    {
        typedef cyPointT Vec3d;

        bool TriangleAABBIntersectTest(const Vec3d& v0, const Vec3d& v1, const Vec3d& v2, const XR::BoundingBox& bbox)
        {
            // 我认为，这里的不等号加上等于号之后，可以作为开集的相交测试
            Vec3d c = bbox.center<Vec3d>();
            Vec3d e = bbox.diagonal<Vec3d>()*0.5;

            Vec3d v00 = v0 - c;
            Vec3d v10 = v1 - c;
            Vec3d v20 = v2 - c;
            //Compute edge vector 
            Vec3d f0 = v10 - v00;
            f0 = Vec3d(fabs(f0[0]), fabs(f0[1]), fabs(f0[2]));
            Vec3d f1 = v20 - v10;
            f1 = Vec3d(fabs(f1[0]), fabs(f1[1]), fabs(f1[2]));
            Vec3d f2 = v00 - v20;
            f2 = Vec3d(fabs(f2[0]), fabs(f2[1]), fabs(f2[2]));
            //Test axes a00 edge-edge test 
            double p0 = v00[2] * v10[1] - v00[1] * v10[2];
            double p2 = v20[2] * (v10[1] - v00[1]) - v20[1] * (v10[2] - v00[2]);
            double r = e[1] * f0[2] + e[2] * f0[1];
            if (std::max(-std::max(p0, p2), std::min(p0, p2)) > r)
                return false;
            // Test axes a01 edge -edge 
            p0 = v10[2] * v20[1] - v10[1] * v20[2];
            p2 = v00[2] * (v20[1] - v10[1]) - v00[1] * (v20[2] - v10[2]);
            r = e[1] * f1[2] + e[2] * f1[1];
            if (std::max(-std::max(p0, p2), std::min(p0, p2)) > r)
                return false;
            // Test axes a02 edge  (dot (v2, a02))
            p0 = v20[2] * v00[1] - v20[1] * v00[2];
            p2 = v10[2] * (v00[1] - v20[1]) - v10[1] * (v00[2] - v20[2]);
            r = e[1] * f2[2] + e[2] * f2[1];
            if (std::max(-std::max(p0, p2), std::min(p0, p2)) > r)
                return false;

            // test axes a10 edge - edge  
            p0 = v00[0] * v10[2] - v00[2] * v10[0];
            p2 = v20[0] * (v10[2] - v00[2]) - v20[2] * (v10[0] - v00[0]);
            r = e[0] * f0[2] + e[2] * f0[0];
            if (std::max(-std::max(p0, p2), std::min(p0, p2)) > r)
                return false;
            p0 = v10[0] * v20[2] - v10[2] * v20[0];
            p2 = v00[0] * (v20[2] - v10[2]) - v00[2] * (v20[0] - v10[0]);
            r = e[0] * f1[2] + e[2] * f1[0];
            if (std::max(-std::max(p0, p2), std::min(p0, p2)) > r)
                return false;
            p0 = v20[0] * v00[2] - v20[2] * v00[0];
            p2 = v10[0] * (v00[2] - v20[2]) - v10[2] * (v00[0] - v20[0]);
            r = e[0] * f2[2] + e[2] * f2[0];
            if (std::max(-std::max(p0, p2), std::min(p0, p2)) > r)
                return false;

            // test axes a20 edge 
            p0 = v00[1] * v10[0] - v00[0] * v10[1];
            p2 = v20[1] * (v10[0] - v00[0]) - v20[0] * (v10[1] - v00[1]);
            r = e[0] * f0[1] + e[1] * f0[0];
            if (std::max(-std::max(p0, p2), std::min(p0, p2)) > r)
                return false;
            p0 = v10[1] * v20[0] - v10[0] * v20[1];
            p2 = v00[1] * (v20[0] - v10[0]) - v00[0] * (v20[1] - v10[1]);
            r = e[0] * f1[1] + e[1] * f1[0];
            if (std::max(-std::max(p0, p2), std::min(p0, p2)) > r)
                return false;
            p0 = v20[1] * v00[0] - v20[0] * v00[1];
            p2 = v10[1] * (v00[0] - v20[0]) - v10[0] * (v00[1] - v20[1]);
            r = e[0] * f2[1] + e[1] * f2[0];
            if (std::max(-std::max(p0, p2), std::min(p0, p2)) > r)
                return false;

            //   /* test in X-direction */
            std::pair<double, double> minmaxValue;
            minmaxValue = std::minmax({ v00[0], v10[0], v20[0] });
            if (minmaxValue.first> e[0] || minmaxValue.second<-e[0]) return false;
            minmaxValue = std::minmax({ v00[1], v10[1], v20[1] });
            if (minmaxValue.first> e[1] || minmaxValue.second<-e[1]) return false;
            minmaxValue = std::minmax({ v00[2], v10[2], v20[2] });
            if (minmaxValue.first> e[2] || minmaxValue.second<-e[2]) return false;

            //test 
            Vec3d normal = (v10 - v00).Cross(v20 - v10);
            double       d = -normal.Dot(v0);
            Vec3d  e1 = e;
            double  r1 = e1.Dot(Vec3d(fabs(normal[0]), fabs(normal[1]), fabs(normal[2])));
            double  s = normal.Dot(c) + d;
            return  (fabs(s) <= (r1));
        }


        inline bool triboxtest(const Triangle* pTri, const XR::BoundingBox& bbox)
        {
            CGALTriangle cgalTri = convertToCGALTriangle(pTri);

            //bool res =  CGAL::do_intersect(CGAL::Bbox_3(bbox.minVal(0), bbox.minVal(1),
            //    bbox.minVal(2), bbox.maxVal(0), bbox.maxVal(1), bbox.maxVal(2)), cgalTri);

            bool res = TriangleAABBIntersectTest(pTri->point(0), pTri->point(1), pTri->point(2), bbox);

            //assert(res == TriangleAABBIntersectTest(pTri->point(0), pTri->point(1), pTri->point(2), bbox));
            return res;
        }

        inline bool fasterTriboxtest(const cyPointT& a, const cyPointT& b, const cyPointT& c, const XR::BoundingBox& bbox)
        {
            double p[3][3] = { { a.x, a.y, a.z }, { b.x, b.y, b.z }, { c.x, c.y, c.z } };
            cyPointT t1 = bbox.center<cyPointT>();
            cyPointT t2 = bbox.diagonal<cyPointT>()*0.5;
            double center[3] = { t1.x, t1.y, t1.z };
            double radius[3] = { t2.x, t2.y, t2.z };
            xrcount++;

            return triBoxOverlap(center, radius, p) == 1;
        }
    }

    int Octree::MAX_TRIANGLE_COUNT = 30;
    int Octree::MAX_LEVEL = 9;

    void Octree::release()
    {
        SAFE_DELETE(mp_root);
        mp_meshes = nullptr;
        m_nMesh = 0;
    }

    Octree::Node* Octree::createRootNode(const Bbox_3& bbox)
    {
        Node* root = new Node;
		root->bbox = XR::BoundingBox(bbox.xmin(), bbox.ymin(), bbox.zmin(),
            bbox.xmax(), bbox.ymax(), bbox.zmax());
        for (uint32_t i = 0; i < m_nMesh; i++)
        {
            auto pcMesh = mp_meshes[i];
			auto &faces = pcMesh->faces();

			for (auto fPtr : faces)
            {
                root->triCount++;
                auto &pTable = root->triTable[i];
                if (!pTable)
                    pTable = new TriList;

				assert(fPtr->degree() == 3);
                pTable->push_back(reinterpret_cast<Triangle*>(fPtr));
            }
        }
        return root;
    }

    void Octree::build(const std::vector<RegularMesh*>& meshList, const Bbox_3& bbox, std::vector<Node*>* isectNodes)
    {
        if (!meshList.size()) throw std::exception("mesh list size 0.");

        mp_meshes = meshList.data();
        m_nMesh = meshList.size();
        mp_root = createRootNode(bbox);

        build(mp_root, 0, isectNodes);
    }

    void Octree::build(Node* root, size_t level, std::vector<Node*>* isectNodes)
    {
        assert(root);

        if (root->triTable.size() <= 1)
        {
            root->type = NODE_SIMPLE;
        }
        else if (root->triCount <= MAX_TRIANGLE_COUNT || level > MAX_LEVEL)
        {
            root->type = NODE_COMPOUND;
            if (isectNodes) isectNodes->push_back(root);
        }
        else
        {
            root->type = NODE_MIDSIDE;
            root->pChildren = new Node[8];

            cyPointT minOffset, maxOffset;
            NodeShape bbox = root->bbox;
            cyPointT step = (bbox.max<cyPointT>() - bbox.min<cyPointT>()) * 0.5;
    
            for (unsigned i = 0; i < 8 ; i++)
            {
                auto pChild = &root->pChildren[i];

                maxOffset = cyPointT(i & 4 ? 0 : -step[0],
                    i & 2 ? 0 : -step[1],
                    i & 1 ? 0 : -step[2]);

                minOffset = cyPointT(i & 4 ? step[0] : 0,
                    i & 2 ? step[1] : 0,
                    i & 1 ? step[2] : 0);

                cyPointT tmpa, tmpb;
                tmpa = bbox.min<cyPointT>() + minOffset;
                tmpb = bbox.max<cyPointT>() + maxOffset;
                pChild->bbox = NodeShape(tmpa.x, tmpa.y, tmpa.z,
                    tmpb.x, tmpb.y, tmpb.z);

                pChild->pParent = root;
            }

            for (auto &triTab: root->triTable)
            {
                uint32_t meshId = triTab.first;
                RegularMesh* pMesh = mp_meshes[meshId];
                TriList &parentMeshes = *triTab.second;

                const size_t tn = parentMeshes.size();
                for (size_t i = 0; i < tn; i++)
                {
                    auto fh = parentMeshes[i];
                    const cyPointT* pts[3] = { &fh->point(0), &fh->point(1), &fh->point(2) };

                    for (unsigned j = 0; j < 8; j++)
                    {
                        Node &child = root->pChildren[j];
                        TriList* triList = nullptr;
                        int pointCount = 0; // -1 no intersection, 1 bbox in box, 0 bbox not in box

                        if (child.bbox.isInclude_left_open_right_close(*pts[0])) pointCount++;
                        if (child.bbox.isInclude_left_open_right_close(*pts[1])) pointCount++;
                        if (child.bbox.isInclude_left_open_right_close(*pts[2])) pointCount++;

                        if (pointCount > 0 || fasterTriboxtest(*pts[0], *pts[1], *pts[2], child.bbox))
                        {
                            if (!triList)
                                triList = child.triTable.emplace(meshId, new TriList()).first->second;

                            triList->push_back(fh);
                            child.triCount++;

                            if (pointCount == 3)
                                break;
                        }
                    }
                }
            }

            root->triTable.clear();
            for (size_t i = 0; i < 8; i++)
                build(root->pChildren + i, level + 1, isectNodes);

        }
    }


} // namespace Boolean

