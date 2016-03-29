#include <vector>

#include "csgdefs.h"
#include "ItstAlg.h"
#include "plane_reps.h"
#include "AlgUserData.h"

namespace
{
    using namespace CSG;

    inline int edge_idx(PosTag tag)
    {
        switch (tag)
        {
        case EDGE_0:
            return 0;
        case EDGE_1:
            return 1;
        case EDGE_2:
            return 2;
        default:
        {
            char _errline[20];
            throw std::exception(ReportErrorString.c_str());
        }
        }
    }

    inline int vertex_idx(PosTag tag)
    {
        switch (tag)
        {
        case VER_0:
            return 0;
        case VER_1:
            return 1;
        case VER_2:
            return 2;
        default:
        {
            char _errline[20];
            throw std::exception(ReportErrorString.c_str());
        }
        }
    }

    inline PosTag edge_tag(int idx)
    {
        switch (idx)
        {
        case 0: return EDGE_0;
        case 1: return EDGE_1;
        case 2: return EDGE_2;
        default: return NONE;
        }
    }

    inline PosTag vertex_tag(int idx)
    {
        switch (idx)
        {
        case 0: return VER_0;
        case 1: return VER_1;
        case 2: return VER_2;
        default: return NONE;
        }
    }

    inline bool is_edge(PosTag tag)
    {
        return tag >= EDGE_0 && tag <= EDGE_2;
    }

    inline bool is_vertex(PosTag tag)
    {
        return tag >= VER_0 && tag <= VER_2;
    }

    inline bool is_same_edge(PosTag tag0, PosTag tag1, int& eId)
    {
        if (tag0 == tag1)
            if (is_edge(tag0))
            {
                eId = edge_idx(tag0);
                return true;
            }

        if (is_edge(tag0) && is_vertex(tag1))
        {
            if (edge_idx(tag0) != vertex_idx(tag1))
            {
                eId = edge_idx(tag0);
                return true;
            }
            else return false;
        }
        else if (is_edge(tag1) && is_vertex(tag0))
        {
            if (edge_idx(tag1) != vertex_idx(tag0))
            {
                eId = edge_idx(tag1);
                return true;
            }
            else return false;
        }
        else if (is_vertex(tag0) && is_vertex(tag1))
        {
            eId = (std::max(vertex_idx(tag0), vertex_idx(tag1)) + 1) % 3;
            return true;
        }
        else
        {
            return false;
        }
    }


    /* A, B 的方向跟cross(ref, triangle)的方向一致*/
    template <class _R>
    static Sign compute_intervals_isectline(CGAL::Oriented_side d[3],
        const CGAL::Triangle_3<_R>& triangle, PBTriangle<_R> &pr, const Plane_ext<_R>& ref,
        PosTag& tagA, PosTag& tagB, Plane_ext<_R>& posA, Plane_ext<_R>& posB)
    {
        std::vector<int> zeroCount;
        std::vector<int> posCount;
        std::vector<int> negCount;
        for (int i = 0; i < 3; i++)
        {
            switch (d[i])
            {
            case CGAL::ON_ORIENTED_BOUNDARY:
                zeroCount.push_back(i);
                break;
            case CGAL::ON_POSITIVE_SIDE:
                posCount.push_back(i);
                break;
            case CGAL::ON_NEGATIVE_SIDE:
                negCount.push_back(i);
                break;
            default:
                ReportError();
                break;
            }
        }

        if (zeroCount.size() == 0)
        {
            int isoVId = (posCount.size() == 1 ? posCount[0] : negCount[0]);
            posA = pr.m_bps[(isoVId + 2) % 3];
            posB = pr.m_bps[(isoVId + 1) % 3];
            tagA = edge_tag((isoVId + 2) % 3);
            tagB = edge_tag((isoVId + 1) % 3);

            if (d[isoVId] == CGAL::ON_NEGATIVE_SIDE)
            {
                std::swap(posA, posB);
                std::swap(tagA, tagB);
            }
            return INTERSECT_ON_LINE;
        }
        else if (zeroCount.size() == 1)
        {
            if (posCount.size() == 0)
                return INTERSECT_ON_POINT;

            int posId = posCount[0];
            int idA = (posId + 2) % 3;
            int idB = (posId + 1) % 3;
            posA = pr.m_bps[idA];
            posB = pr.m_bps[idB];

            if (idA == zeroCount[0])
            {
                tagA = edge_tag(idA);
                tagB = vertex_tag(idA);
            }
            else
            {
                tagB = edge_tag(idA);
                tagA = vertex_tag(idA);
            }
            return INTERSECT_ON_LINE;
        }
        else if (zeroCount.size() == 2)
        {
            int isoVId = (posCount.size() == 1 ? posCount[0] : negCount[0]);
            posA = pr.m_bps[(isoVId + 2) % 3];
            posB = pr.m_bps[(isoVId + 1) % 3];

            tagA = vertex_tag((isoVId + 1) % 3);
            tagB = vertex_tag((isoVId + 2) % 3);

            if (d[isoVId] == CGAL::ON_NEGATIVE_SIDE)
            {
                std::swap(posA, posB);
                std::swap(tagA, tagB);
            }
            return INTERSECT_ON_LINE;
        }
        else
        {
            ReportError();
            return NOT_INTERSECT;
        }
    }

    template <class _R>
    void makePositive(const CGAL::Plane_3<_R>& p, const CGAL::Plane_3<_R>& q, CGAL::Plane_3<_R>& input)
    {
        if (CGAL::determinant(p.orthogonal_vector(),
            q.orthogonal_vector(), input.orthogonal_vector()) < 0.f)
            input = input.opposite();
    }

    template <class _R>
    Sign tri_tri_intersect(const CGAL::Triangle_3<_R> t[2], const Plane_ext<_R> p[2],
        TriTriIsectResult<_R>* result, MyMesh::Face_handle fhs[2])
    {
        typedef _R::FT FT;
        typedef _R::Vector_3 Vec3d;

        CGAL::Oriented_side db[3];
        for (size_t i = 0; i < 3; i++)
            db[i] = p[0].oriented_side(t[1].vertex(i));

        if ((db[0] == db[1]) && (db[1] == db[2]))
        {
            if (db[0] == CGAL::ON_ORIENTED_BOUNDARY) return COPLANAR;
            else return NOT_INTERSECT;
        }

        CGAL::Oriented_side da[3];
        for (size_t i = 0; i < 3; i++)
            da[i] = p[1].oriented_side(t[0].vertex(i));

        if ((da[0] == da[1]) && (da[1] == da[2]))
        {
            if (da[0] == CGAL::ON_ORIENTED_BOUNDARY) return COPLANAR;
            else return NOT_INTERSECT;
        }

        for (int i = 0; i < 2; i++)
        {
            if (!fhs[i]->data->planeRep)
                fhs[i]->data->planeRep = new PBTriangle<_R>(t[i]);
        }

        PosTag tagA[2], tagB[2];
        Plane_ext<_R> posA[2], posB[2];
        Sign sign;

        // 规定cross(n1, n0)为正方向
        sign = compute_intervals_isectline(da, t[0], *fhs[0]->data->planeRep, p[1], tagA[0], tagB[0], posA[0], posB[0]);
        if (sign == INTERSECT_ON_POINT)
            return sign;

        // 规定cross(n0, n1)为正方向
        sign = compute_intervals_isectline(db, t[1], *fhs[1]->data->planeRep, p[0], tagA[1], tagB[1], posA[1], posB[1]);
        if (sign == INTERSECT_ON_POINT)
            return sign;

        // 统一正方向cross(n0, n1)
        std::swap(tagA[0], tagB[0]);
        std::swap(posA[0], posB[0]);

        makePositive(p[0], p[1], posA[0]);
        makePositive(p[0], p[1], posA[1]);
        makePositive(p[0], p[1], posB[0]);
        makePositive(p[0], p[1], posB[1]);

        assert(orientation(p[0], p[1], posB[0], posA[0]) > 0);
        assert(orientation(p[0], p[1], posB[1], posA[1]) > 0);

        double cmpA0B1 = orientation(p[0], p[1], posB[1], posA[0]);
        if (cmpA0B1 < 0.) return NOT_INTERSECT;
        if (cmpA0B1 == 0.) return INTERSECT_ON_POINT;

        double cmpA1B0 = orientation(p[0], p[1], posB[0], posA[1]);
        if (cmpA1B0 < 0.) return NOT_INTERSECT;
        if (cmpA1B0 == 0.) return INTERSECT_ON_POINT;

        double Acmp = orientation(p[0], p[1], posA[0], posA[1]);
        double Bcmp = orientation(p[0], p[1], posB[0], posB[1]);

        //static int count = 0;
        //std::cout << count++ << std::endl;
        //std::cout << t[0] << std::endl;
        //std::cout << t[1] << std::endl;

        // if (Acmp > 0): A0在A1右边，取A0
        result->tagA[0] = Acmp >= 0.0 ? tagA[0] : INNER;
        result->tagA[1] = Acmp <= 0.0 ? tagA[1] : INNER;
        result->A = PBPoint<_R>(p[0], p[1], (Acmp >= 0.0 ? posA[0] : posA[1]));

        // if (Bcmp > 0): B0在B1右边，取B1
        result->tagB[0] = Bcmp <= 0.0 ? tagB[0] : INNER;
        result->tagB[1] = Bcmp >= 0.0 ? tagB[1] : INNER;
        result->B = PBPoint<_R>(p[0], p[1], (Bcmp >= 0.0 ? posB[1] : posB[0]));

        //result->A.computeCoord();
        //result->B.computeCoord();
        //std::cout << result->A.coord[0] << " " << result->A.coord[1] << " " << result->A.coord[2] << " " << std::endl;
        //std::cout << result->B.coord[0] << " " << result->B.coord[1] << " " << result->B.coord[2] << " " << std::endl;

        assert(result->tagA[0] != result->tagB[0]);
        assert(result->tagA[1] != result->tagB[1]);

        return INTERSECT_ON_LINE;
    }
}


namespace CSG
{
    void AdjacentGraph::getIntersectPrimitives(int meshId, std::vector<int>& prims)
    {
        for (size_t i = 0; i < m_sz; i++)
        {
            if (getValue(meshId, i))
                prims.push_back(i);
        }
        return;
    }

    ItstAlg::ItstAlg(std::vector<MyMesh*>* meshes)
    {
        mp_meshList = meshes;
    }


    ItstAlg::~ItstAlg()
    {
        SAFE_DELETE(mp_adjGraph);
    }


    void ItstAlg::doIntersection(std::vector<Octree::Node*>& intersectLeaves)
    {
        MeshIdTriIdMap antiOverlapMap;
        antiOverlapMap.max_load_factor(0.6f);
        mp_adjGraph = new AdjacentGraph(mp_meshList->size());
        auto &meshList = *mp_meshList;

        for (Octree::Node* leaf : intersectLeaves)
        {
            auto iEnd = leaf->triTable.cend();
            decltype(leaf->triTable.begin()) meshItr[2];
            for (meshItr[0] = leaf->triTable.begin(); meshItr[0] != iEnd; ++meshItr[0])
            {
                meshItr[1] = meshItr[0]; ++meshItr[1];
                for (; meshItr[1] != iEnd; ++meshItr[1])
                {
                    uint32_t meshId[2] = { meshItr[0]->first, meshItr[1]->first };
                    MyMesh* meshes[2] = { meshList[meshId[0]], meshList[meshId[1]] };

                    // 这里假设map的遍历是保序的，因此meshIdPair自动的分为大小
                    IndexPair meshIdPair;
                    MakeIndex(meshId, meshIdPair);

                    TriIdSet* antiOverlapSet = nullptr;
                    auto searchRes = antiOverlapMap.find(meshIdPair);
                    if (searchRes == antiOverlapMap.end())
                    {
                        antiOverlapSet = new TriIdSet;
                        antiOverlapSet->max_load_factor(0.6f);
                        antiOverlapMap.emplace(meshIdPair, antiOverlapSet);
                    }
                    else antiOverlapSet = searchRes->second;

                    for (MyMesh::Face_handle fh0 : *meshItr[0]->second)
                    {
                        for (MyMesh::Face_handle fh1 : *meshItr[1]->second)
                        {
                            uint32_t triId[2] = { fh0->id(), fh1->id() };
                            IndexPair triIdPair;
                            MakeIndex(triId, triIdPair);

                            if (antiOverlapSet->find(triIdPair) != antiOverlapSet->end())
                                continue;

                            antiOverlapSet->insert(triIdPair);

                            if (IntersectionTest(fh0, fh1, antiOverlapSet, meshId))
                                mp_adjGraph->setValue(meshId[0], meshId[1], true);
                        }
                    }
                }
            }
        }
    }

    bool ItstAlg::checkManifoldEdge(FH fh0, FH fh1, TriIdSet* overlaps, TriTriIsectResult<K> &result, int res[])
    {
        K::Point_3 vthiz, vthat;
        K::Plane_3 plane;
        FH fh[2] = { fh0, fh1 };
        res[0] = false; res[1] = false;
        bool detected = false;
        
        for (int i = 0; i < 2; i++)
        {
            // 首先检查存不存在交点都在同一个边界上的情况
            int edgeId = -1;
            auto &out = res[i];
            out = 2; // add lines
            if (is_same_edge(result.tagA[i], result.tagB[i], edgeId))
            {
                out = 0; // add vertices
                detected = true;
                int idx = edgeId;
                vthiz = fh[i]->data->triangle.vertex(idx);
                MyMesh::Halfedge_handle oppoHE = fh[i]->edges[idx]->opposite();
                vthat = oppoHE->facet()->data->triangle.vertex(oppoHE->id);
                plane = fh[(i+1)%2]->data->planeRep->get_sp();

                uint32_t triId[2];
                triId[i] = oppoHE->facet()->id();
                triId[(i + 1) % 2] = fh[(i + 1) % 2]->id();
                IndexPair triIdPair;
                MakeIndex(triId, triIdPair);
                overlaps->insert(triIdPair);

                if ((plane.has_on_positive_side(vthiz) && plane.has_on_positive_side(vthat)) ||
                    (plane.has_on_negative_side(vthiz) && plane.has_on_negative_side(vthat)))
                    out = 0; // add none
            }
        }

        return detected;
    }

    int ItstAlg::checkDuplicatedPoints(std::vector<VProxyItr>& plist, PBPoint<K> point, VProxyItr& proxy)
    {
        for (size_t i = 0; i < plist.size(); i++)
        {
            if (plist[i].pointer()->pos == point)
            {
                proxy = plist[i];
                return i;
            }
        }
        return -1;
    }

    int ItstAlg::checkDuplicatedPoints(PBPoint<K> point, FH fhs, PosTag tags, VProxyItr& outcome)
    {
        int oId = -1;
        switch (tags)
        {
        case INNER:
            oId = checkDuplicatedPoints(fhs->data->itstTri->inVertices, point, outcome);
            break;
        case VER_0:
        case VER_1:
        case VER_2:
            if (fhs->vertices[vertex_idx(tags)]->data)
                oId = 0;
            break;
        case EDGE_0:
        case EDGE_1:
        case EDGE_2:
        {
            auto &pEdge = fhs->edges[edge_idx(tags)]->data;
            if (!pEdge)
            {
                pEdge.reset(new UserEData);
                fhs->edges[edge_idx(tags)]->opposite()->data = pEdge;
            }
            else oId = checkDuplicatedPoints(pEdge->vertices, point, outcome);
            break;
        }
        default:
            ReportError();
            break;
        }
        return oId;
    }

    int ItstAlg::addPoint(FH fh, PosTag tags, VProxyItr proxy)
    {
        switch (tags)
        {
        case INNER:
            fh->data->itstTri->inVertices.push_back(proxy);
            return fh->data->itstTri->inVertices.size() - 1;
        case VER_0:
        case VER_1:
        case VER_2:
            assert(!fh->vertices[vertex_idx(tags)]->data);
            fh->vertices[vertex_idx(tags)]->data.reset(new UserVData(&proxy));
            return 0;
        case EDGE_0:
        case EDGE_1:
        case EDGE_2:
            fh->edges[edge_idx(tags)]->data->vertices.push_back(proxy);
            return fh->edges[edge_idx(tags)]->data->vertices.size() - 1;
        default:
            ReportError();
            break;
        }
    }

    void ItstAlg::getVProxy(PBPoint<K> point, int addwhat[2], FH fhs[2], PosTag tags[2], 
        int oId[2], VProxyItr outproxy[2], uint32_t meshId[2])
    {
        VProxyItr proxy;

        for (int i = 0; i < 2; i++)
        {
            if (addwhat[i] > 0)
                oId[i] = checkDuplicatedPoints(point, fhs[i], tags[i], outproxy[i]);
        }

        if (oId[0] == -1 && oId[1] == -1)
        {
            vEnt.emplace_back(new VEntity);
            auto &entity = vEnt.back();
            entity->pos = point;
            for (int i = 0; i < 2; i++)
            {
                if (addwhat[i] > 0)
                    entity->addContext(meshId[i], fhs[i], tags[i]);
            }

            auto veItr = vEnt.end(); veItr--;
            vProxy.push_back(veItr);
            proxy = vProxy.end(); proxy--;
            outproxy[0] = outproxy[1] = proxy;
        }
        else if (oId[0] > -1 && oId[1] > -1)
        {
            if (*outproxy[0].pointer() < *outproxy[1].pointer())
                mergeProxy(outproxy[0], outproxy[1]);
            else
                mergeProxy(outproxy[1], outproxy[0]);
        }
        else if (oId[0] == -1)
            outproxy[0] = outproxy[1];
        else
            outproxy[1] = outproxy[0];
    }

    bool ItstAlg::IntersectionTest(FH fh0, FH fh1, TriIdSet* overlaps, uint32_t meshId[2])
    {
        K::Triangle_3 t[2] = { fh0->data->triangle, fh1->data->triangle };
        Plane_ext<K> sp[2] = { fh0->data->sp, fh1->data->sp };
        FH fhs[2] = { fh0, fh1 };
        TriTriIsectResult<K> result;

        Sign sign = tri_tri_intersect(t, sp, &result, fhs);

        if (sign == NOT_INTERSECT || sign == INTERSECT_ON_POINT || sign == COPLANAR)
            return false;

        int addwhat[2]; // 0: nothing, 1: add vertex, 2: add line
        checkManifoldEdge(fh0, fh1, overlaps, result, addwhat);

        if (addwhat[0] == 0 && addwhat[1] == 0)
            return false;

        ItstTriangle*& it0 = fh0->data->itstTri;
        ItstTriangle*& it1 = fh1->data->itstTri;

        if (!it0) it0 = new ItstTriangle(fh0);
        if (!it1) it1 = new ItstTriangle(fh1);

        int oIdA[2] = { -1, -1 }, oIdB[2] = { -1, -1 };
        VProxyItr proxyA[2], proxyB[2];

        getVProxy(result.A, addwhat, fhs, result.tagA, oIdA, proxyA, meshId);
        getVProxy(result.B, addwhat, fhs, result.tagB, oIdB, proxyB, meshId);
        
        for (int i = 0; i < 2; i++)
        {
            if (addwhat[i] > 0)
            {
                if (oIdA[i] == -1)
                    oIdA[i] = addPoint(fhs[i], result.tagA[i], proxyA[i]);

                if (oIdB[i] == -1)
                    oIdB[i] = addPoint(fhs[i], result.tagB[i], proxyB[i]);
            }
        }

        // prepare the line segments
        ItstLine line;
        for (int i = 0; i < 2; i++)
        {
            if (addwhat[i] > 1)
            {
                line.pts[0].tag = result.tagA[i];
                line.pts[1].tag = result.tagB[i];

                line.pts[0].idx = oIdA[i];
                line.pts[1].idx = oIdB[i];

                fhs[i]->data->itstTri->isectLines.push_back(line);
                fhs[i]->data->itstTri->meshIds.insert(meshId[(i + 1) % 2]);
            }
        }
        return true;
    }


    void ItstAlg::mergeProxy(VProxyItr a, VProxyItr b)
    {
        if (a.pointer() == b.pointer()) return;

        a.pointer()->ctx.insert(a.pointer()->ctx.end(), 
            b.pointer()->ctx.begin(), b.pointer()->ctx.end());
        *b = *a;
    }
}