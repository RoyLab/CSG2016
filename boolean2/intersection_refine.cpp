#include "precompile.h"
#include "intersection.h"
#include "hybrid_geometry.h"

namespace Boolean
{
    void EdgeInsctData::refine(void* pData)
    {
        if (isRefined()) return;

        std::vector<PlaneVertex> seqs(points.size());
        auto vItr = points.begin();
        for (int i = 0; i < points.size(); i++, vItr++)
            seqs[i].id = *vItr;

        std::map<VertexIndex, XPlane> v2p;
        for (auto& set : inscts)
        {
            for (auto &pbi : set.second)
            {
                v2p[pbi.ends[0]] = pbi.pends[0];
                v2p[pbi.ends[1]] = pbi.pends[1];
            }
        }

        EdgeAuxiliaryStructure data = *(EdgeAuxiliaryStructure*)(pData);
        // assign each vertex a plane and the corresponding id
        for (int i = 0; i < points.size(); i++)
        {
            auto res0 = v2p.find(seqs[i].id);
            if (res0 == v2p.end())
            {
                // some vertex cannot be assigned a valid plane
                // therefore we should find it manually from Plane Triples
                auto &vRef = xvertex(seqs[i].id);
                if (vRef.isPlaneRep())
                {
                    auto& xpointRef = vRef.ppoint();
                    for (int j = 0; j < 3; j++)
                    {
                        Real fres = data.line.dot(xpointRef.plane(j));
                        if (fres == Real(0)) continue;

                        if (fres > 0)
                            seqs[i].plane = xpointRef.plane(j);
                        else if (fres < 0)
                            seqs[i].plane = xpointRef.plane(j).opposite();
                        break;
                    }
                }
            }
            else seqs[i].plane = res0->second;
        }

        LinOrderObj orderObj = { data.line };
        std::sort(seqs.begin(), seqs.end(), orderObj);
        std::vector<EdgePBI> newPbi(points.size() + 1);
        std::map<VertexIndex, uint32_t> idmap;
        idmap[data.start] = 0;
        idmap[data.end] = points.size() + 1;
        for (int i = 0; i < seqs.size(); i++)
        {
            idmap[seqs[i].id] = i + 1;
            newPbi[i].ends[1] = seqs[i].id;
            newPbi[i + 1].ends[0] = seqs[i].id;
            newPbi[i].pends[1] = seqs[i].plane;
            newPbi[i + 1].pends[0] = seqs[i].plane;
        }
        newPbi[0].ends[0] = data.start;
        newPbi[points.size()].ends[1] = data.end;

        // add pbi intersections' neighborInfo into newly constructed pbi.
        for (auto& set : inscts)
        {
            for (auto &pbi : set.second)
            {
                assert(idmap.find(pbi.ends[0]) != idmap.end());
                assert(idmap.find(pbi.ends[1]) != idmap.end());

                uint32_t start = idmap[pbi.ends[0]];
                uint32_t last = idmap[pbi.ends[1]];

                assert(start < last);
                for (int i = start; i < last; i++)
                {
                    assert(pbi.neighbor.size() == 1);
                    newPbi[i].neighbor.push_back(*pbi.neighbor.begin());
                }
            }
        }
        inscts.clear();
        auto& slot = inscts[INVALID_UINT32];
        for (int i = 0; i < newPbi.size(); i++)
            slot.push_back(newPbi[i]);

        bRefined = true;
    }

    // strayVertices是孤立点，从其他的三面交点处传过来的？
    bool checkIsolatedPart(Triangle* pTri, std::vector<VertexIndex>& strayVertices)
    {
        struct ColorVertex
        {
            std::set<VertexIndex> neighbors;
            int color = -1;
        } defaultCV;

        FaceInsctData* thiz = pTri->inscts;
        // initialize the graph
        std::map<VertexIndex, ColorVertex> data;
        for (auto& v : thiz->points)
            data[v.vId] = defaultCV;

        for (auto& pbiSet : thiz->inscts)
        {
            for (FacePBI& fPbi : pbiSet.second)
            {
                data[fPbi.ends[0]].neighbors.insert(fPbi.ends[1]);
                data[fPbi.ends[1]].neighbors.insert(fPbi.ends[0]);
            }
        }

        // init vertex set
        std::set<VertexIndex> edgeVertices;
        for (int i = 0; i < 3; i++)
        {
            edgeVertices.insert(pTri->vertexId(i));
            MyEdge& edge = pTri->edge(i);
            if (!edge.inscts) continue;
            for (VertexIndex vId : edge.inscts->points)
                edgeVertices.insert(vId);
        }

        // colorization
        int color = 0;
        uint32_t isoCount = MAX_MESH_COUNT;
        bool isIsolated = false, boolRes = false;
        std::stack<VertexIndex> vStack;
        std::vector<VertexIndex> history;
        for (auto& vItem : data)
        {
            if (vItem.second.color == -1)
            {
                ++color;
                isIsolated = true;
                assert(vStack.empty());
                vStack.push(vItem.first);
                while (!vStack.empty())
                {
                    VertexIndex curVId = vStack.top();
                    vStack.pop();
                    auto& curV = data[curVId];
                    if (curV.color != -1) continue;

                    history.push_back(curVId);
                    if (edgeVertices.find(curVId) != edgeVertices.end())
                        isIsolated = false;

                    curV.color = color;
                    for (VertexIndex vId : curV.neighbors)
                    {
                        if (data[vId].color == -1)
                            vStack.push(vId);
                    }
                }

                if (isIsolated)
                {
                    VertexIndex chooseVertex = INVALID_UINT32;
                    MyEdge::SIndex eId = -1;
                    for (VertexIndex vId : history)
                    {
                        for (auto &v : thiz->points)
                        {
                            if (v.vId == vId)
                            {
                                chooseVertex = vId;
                                eId = v.eId;
                                break;
                            }
                        }

                        if (eId >= 0) break;
                    }
                    assert(eId != -1);
                    if (eId >= 0)
                    {
                        FacePBI fPbi;
                        fPbi.ends[0] = pTri->vertexId(0);
                        fPbi.ends[1] = chooseVertex;

                        MyEdge& crossEdgeRef = xedge(eId);
                        fPbi.vertPlane = XPlane(pTri->vertex(0).point(),
                            xvertex(crossEdgeRef.ends[0]).point(), xvertex(crossEdgeRef.ends[1]).point());

                        PlaneLine line(pTri->supportingPlane(), fPbi.vertPlane);
                        fPbi.pends[0] = pickPositiveVertical(line, xvertex(fPbi.ends[0]));
                        fPbi.pends[1] = pickPositiveVertical(line, xvertex(fPbi.ends[1]));
                        assert(line.dot(fPbi.pends[0]) > 0);
                        assert(line.dot(fPbi.pends[1]) > 0);
                        if (line.linear_order_unsafe(fPbi.pends[0], fPbi.pends[1]) < 0)
                        {
                            std::swap(fPbi.ends[0], fPbi.ends[1]);
                            std::swap(fPbi.pends[0], fPbi.pends[1]);
                        }

                        assert(thiz->inscts.find(isoCount) == thiz->inscts.end());
                        thiz->inscts[isoCount++].push_back(fPbi);
                    }
                    else
                    {
                        assert(eId == -2 && history.size() == 1);
                        strayVertices.push_back(history[0]);
                    }
                    boolRes = true;
                }
                history.clear();
            }
        }

        return boolRes;
    }

    void removeOverlapPBI(FaceInsctData * thiz)
    {
        std::map<VertexIndex, std::set<FacePBI*>> data;
        std::vector<decltype(thiz->inscts[0].begin())> garbage;
        for (auto &pbiSet : thiz->inscts)
        {
            for (auto pbiItr = pbiSet.second.begin(); pbiItr != pbiSet.second.end(); ++pbiItr)
            {
                FacePBI& pbi = *pbiItr;
                bool found = false;
                for (auto alreadyHere : data[pbi.ends[0]]) // search in current graph
                {
                    if (alreadyHere->ends[0] == pbi.ends[1] ||
                        alreadyHere->ends[1] == pbi.ends[1]) // if has
                    {
                        for (NeighborInfo& nInfo : pbi.neighbor)
                        {
                            bool unique = true; // search if the neighInfo is already there, == std::find
                            for (NeighborInfo& nInfo2 : alreadyHere->neighbor)
                            {
                                if (nInfo.neighborMeshId == nInfo.neighborMeshId)
                                {
                                    unique = false;
                                    break;
                                }
                            }

                            if (unique)
                            {
                                alreadyHere->neighbor.push_back(nInfo);
                            }
                        }
                        found = true;
                        break;
                    }
                }

                if (found)
                {
                    garbage.push_back(pbiItr);
                    continue;
                }

                data[pbi.ends[0]].insert(&pbi);
                data[pbi.ends[1]].insert(&pbi);
            }

            for (auto& pbiItr : garbage)
                pbiSet.second.erase(pbiItr);
            garbage.clear();
        }
    }

    IndexPair makePbiIndex(const PBIRep* rep)
    {
        IndexPair res;
        MakeIndex(rep->ends, res);
        return res;
    }

    void FaceInsctData::resolveIntersection(Triangle* pTri, std::vector<VertexIndex>* strayVertices)
    {
        XPlane triSp = pTri->supportingPlane();
        std::map<IndexPair, std::shared_ptr<FacePBITessData>> tessData;
        for (auto setItr = inscts.begin(); setItr != inscts.end(); ++setItr)
        {
            auto setItr2 = setItr; ++setItr2;
            for (; setItr2 != inscts.end(); ++setItr2)
            {
                for (auto pbiItr = setItr->second.begin(); pbiItr != setItr->second.end(); ++pbiItr)
                {
                    FacePBI& pbi = *pbiItr;
                    for (auto pbiItr2 = setItr2->second.begin(); pbiItr2 != setItr2->second.end(); ++pbiItr2)
                    {
                        PlaneVertex slots[2][2]; // 需要添加的顶点
                        FacePBI& pbi2 = *pbiItr2;
                        Oriented_side side[2][2];
                        side[0][0] = orientation(pbi2.vertPlane, pbi.ends[0]);
                        side[0][1] = orientation(pbi2.vertPlane, pbi.ends[1]);
                        if (side[0][0] == side[0][1])
                        {
                            if (side[0][0] == ON_ORIENTED_BOUNDARY) // 共线情况
                            {
                                PlaneLine line(triSp, pbi.vertPlane);
                                assert(line.dot(pbi2.vertPlane) == 0.);
                                Real dotRes = line.dot(pbi2.pends[0]);
                                assert(dotRes != 0.);
                                bool inverseLine = dotRes > 0 ? false : true;

                                PlaneVertex pbi2ends[2] = { { pbi2.pends[0], pbi2.ends[0] },{ pbi2.pends[1], pbi2.ends[1] } };
                                if (inverseLine)
                                {
                                    pbi2ends[0].plane.inverse();
                                    pbi2ends[1].plane.inverse();
                                    std::swap(pbi2ends[0], pbi2ends[1]);
                                }

                                // linear order 计算overlap
                                assert(line.linear_order(pbi2ends[0].plane, pbi2ends[1].plane) > 0);
                                assert(line.linear_order_unsafe(pbi2ends[0].plane, pbi2ends[1].plane) > 0);
                                if (line.linear_order_unsafe(pbi2ends[0].plane, pbi.pends[1]) <= 0 ||
                                    line.linear_order_unsafe(pbi.pends[0], pbi2ends[1].plane) <= 0)
                                    continue;

                                int compRes[2] = {
                                    line.linear_order_unsafe(pbi.pends[0], pbi2ends[0].plane),
                                    line.linear_order_unsafe(pbi.pends[1], pbi2ends[1].plane),
                                };

                                assert(!(compRes[0] == 0 && compRes[1] == 0));

                                if (compRes[0] > 0)
                                    slots[0][0] = pbi2ends[0];
                                else if (compRes[0] < 0)
                                {
                                    slots[1][0].id = pbi.ends[0];
                                    slots[1][0].plane = pbi.pends[0];
                                    if (inverseLine) slots[1][0].plane.inverse();
                                }

                                if (compRes[1] > 0)
                                {
                                    slots[1][1].id = pbi.ends[1];
                                    slots[1][1].plane = pbi.pends[1];
                                    if (inverseLine) slots[1][1].plane.inverse();
                                }
                                else if (compRes[1] < 0)
                                    slots[0][1] = pbi2ends[1];
                            }
                            else continue;
                        }
                        else
                        {
                            side[1][0] = orientation(pbi.vertPlane, pbi2.ends[0]);
                            side[1][1] = orientation(pbi.vertPlane, pbi2.ends[1]);

                            if (side[1][0] == side[1][1])
                            {
                                assert(side[1][0] != ON_ORIENTED_BOUNDARY);//如果是共线情况，在前面应当已经被测试到
                                continue; // 不相交
                            }

                            if (side[0][0] * side[0][1] == 0)
                            {
                                if (side[1][0] * side[1][1] == 0) // 相交于某个顶点，不用split
                                    continue;

                                int addedTarget = side[0][0] == ON_ORIENTED_BOUNDARY ? 0 : 1;
                                slots[1][0].id = pbi.ends[addedTarget];
                                slots[1][0].plane = pbi.vertPlane;
                                PlaneLine(triSp, pbi2.vertPlane).make_positive(slots[1][0].plane);
                            }
                            else
                            {
                                assert(side[0][0] * side[0][1] == -1);
                                if (side[1][0] * side[1][1] == 0)
                                {
                                    int addedTarget = side[1][0] == ON_ORIENTED_BOUNDARY ? 0 : 1;
                                    slots[0][0].id = pbi2.ends[addedTarget];
                                    slots[0][0].plane = pbi2.vertPlane;
                                    PlaneLine(triSp, pbi.vertPlane).make_positive(slots[0][0].plane);
                                }
                                else
                                {
                                    assert(side[1][0] * side[1][1] == -1);
                                    // new vertex
                                    XPlane thirdPlane = pbi2.vertPlane;
                                    PlaneLine(triSp, pbi.vertPlane).make_positive(thirdPlane); // 似乎不需要，可以尝试注释这一句
                                    PlanePoint newPoint(triSp, pbi.vertPlane, thirdPlane);

                                    uint32_t* newPos;
                                    newPos = point(newPoint, -1);
                                    std::vector<uint32_t*> vecs;
                                    vecs.push_back(newPos);
                                    uint32_t minVal = *newPos;

                                    // 去所有的邻居看一看
                                    FacePBI* twoPbiPtr[2] = { &pbi, &pbi2 };

                                    for (int i = 0; i < 2; i++)
                                    {
                                        int i2 = (i + 1) % 2;
                                        for (NeighborInfo& neiInfo : twoPbiPtr[i]->neighbor)
                                        {
                                            if (neiInfo.type == NeighborInfo::Edge)
                                            {
                                                MyEdge& cur_con = xedge(neiInfo.neighborEdgeId);
                                                assert(cur_con.inscts);
                                                newPos = cur_con.inscts->point(newPoint);
                                            }
                                            else
                                            {
                                                assert(neiInfo.type == NeighborInfo::Face);
                                                assert(neiInfo.pTrangle->inscts);
                                                int eId = -1;
                                                if (twoPbiPtr[i2]->neighbor.empty())
                                                    eId = -2;

                                                newPos = neiInfo.pTrangle->inscts->point(newPoint, eId);
                                            }
                                            if (*newPos < minVal)
                                                minVal = *newPos;
                                            vecs.push_back(newPos);
                                        }
                                    }

                                    // 所有的邻居都没有，那就真没有了
                                    if (minVal == INVALID_UINT32)
                                    {
                                        minVal = GlobalData::getObject()->insertVertex(newPoint);
                                    }
                                    for (uint32_t *pInt : vecs)
                                        *pInt = minVal;

                                    slots[0][0].id = minVal;
                                    slots[0][0].plane = pbi2.vertPlane;
                                    PlaneLine(triSp, pbi.vertPlane).make_positive(slots[0][0].plane);

                                    slots[1][0].id = minVal;
                                    slots[1][0].plane = pbi.vertPlane;
                                    PlaneLine(triSp, pbi2.vertPlane).make_positive(slots[1][0].plane);
                                }
                            }

                        }

                        if (slots[0][0].plane.is_valid() || slots[0][1].plane.is_valid())
                        {
                            IndexPair pbiIndex = makePbiIndex(&pbi);
                            auto pbiData = tessData.find(pbiIndex);
                            if (pbiData == tessData.end())
                            {
                                FacePBITessData *tessItem = new FacePBITessData;
                                tessItem->pContainer = &setItr->second;
                                tessItem->ptr = pbiItr;
                                auto insertRes = tessData.insert(decltype(tessData)::value_type(pbiIndex,
                                    std::shared_ptr<FacePBITessData>(tessItem)));

                                assert(insertRes.second);
                                pbiData = insertRes.first;
                            }

                            if (slots[0][0].plane.is_valid())
                                pbiData->second->points.push_back(slots[0][0]);

                            if (slots[0][1].plane.is_valid())
                                pbiData->second->points.push_back(slots[0][1]);
                        }

                        if (slots[1][0].plane.is_valid() || slots[1][1].plane.is_valid())
                        {
                            IndexPair pbiIndex2 = makePbiIndex(&pbi2);
                            auto pbiData = tessData.find(pbiIndex2);
                            if (pbiData == tessData.end())
                            {
                                FacePBITessData *tessItem = new FacePBITessData;
                                tessItem->pContainer = &setItr2->second;
                                tessItem->ptr = pbiItr2;
                                auto insertRes = tessData.insert(decltype(tessData)::value_type(pbiIndex2,
                                    std::shared_ptr<FacePBITessData>(tessItem)));

                                assert(insertRes.second);
                                pbiData = insertRes.first;
                            }

                            if (slots[1][0].plane.is_valid())
                                pbiData->second->points.push_back(slots[1][0]);

                            if (slots[1][1].plane.is_valid())
                                pbiData->second->points.push_back(slots[1][1]);
                        }
                    }
                }
            }
        }

        if (strayVertices)
        {
            for (VertexIndex strayV : *strayVertices)
            {
                if (!xvertex(strayV).isPlaneRep()) continue; // 不合理，但是先这样吧
                for (auto setItr = inscts.begin(); setItr != inscts.end(); ++setItr)
                {
                    for (auto pbiItr = setItr->second.begin(); pbiItr != setItr->second.end(); ++pbiItr)
                    {
                        FacePBI& fPbi = *pbiItr;
                        if (fPbi.vertPlane.has_on(xvertex(strayV).ppoint()))
                        {
                            IndexPair pbiIndex = makePbiIndex(&fPbi);
                            auto pbiData = tessData.find(pbiIndex);
                            if (pbiData == tessData.end())
                            {
                                FacePBITessData *tessItem = new FacePBITessData;
                                tessItem->pContainer = &setItr->second;
                                tessItem->ptr = pbiItr;
                                auto insertRes = tessData.insert(decltype(tessData)::value_type(pbiIndex,
                                    std::shared_ptr<FacePBITessData>(tessItem)));

                                assert(insertRes.second);
                                pbiData = insertRes.first;
                            }

                            PlaneLine line(pTri->supportingPlane(), fPbi.vertPlane);
                            line.pick_positive_vertical_plane(xvertex(strayV).ppoint());
                            pbiData->second->points.push_back(
                                PlaneVertex{ line.pickPositiveVertical(xvertex(strayV).ppoint()), strayV });
                            break; // 如果需要加两个以上的fpbi的话，这个点应该已经被探测出来了
                        }
                    }
                }
            }
        }

        for (auto &pPair : tessData)
        {
            FacePBITessData* pData = pPair.second.get();
            PlaneLine line(triSp, pData->ptr->vertPlane);
            LinOrderObj orderObj = { line };

            auto &inserted = pData->points;
            std::sort(inserted.begin(), inserted.end(), orderObj);
            inserted.erase(std::unique(inserted.begin(), inserted.end(),
                [](const PlaneVertex &a, const PlaneVertex &b)->bool {
                return a.id == b.id;
            }), inserted.end());

            std::vector<FacePBI> newPbi(inserted.size() + 1, *pData->ptr);
            std::map<VertexIndex, uint32_t> idmap;
            idmap[pData->ptr->ends[0]] = 0;
            idmap[pData->ptr->ends[1]] = inserted.size() + 1;
            for (int i = 0; i < inserted.size(); i++)
            {
                idmap[inserted[i].id] = i + 1;
                newPbi[i].ends[1] = inserted[i].id;
                newPbi[i + 1].ends[0] = inserted[i].id;
                newPbi[i].pends[1] = inserted[i].plane;
                newPbi[i + 1].pends[0] = inserted[i].plane;
            }
            newPbi[0].ends[0] = pData->ptr->ends[0];
            newPbi[inserted.size()].ends[1] = pData->ptr->ends[1];

            pData->pContainer->erase(pData->ptr);
            pData->pContainer->insert(pData->pContainer->end(), newPbi.begin(), newPbi.end());
            newPbi.clear();
        }

        removeOverlapPBI(this);
    }


    void FaceInsctData::refine(void* pData)
    {
        removeOverlapPBI(this);
        Triangle* pTri = reinterpret_cast<Triangle*>(pData);
        if (!isRefined() && inscts.size() >= 2)
            resolveIntersection(pTri);

        std::vector<VertexIndex> strayVertices;
        if (checkIsolatedPart(pTri, strayVertices))
            resolveIntersection(pTri, &strayVertices);

        bRefined = true;
    }
}