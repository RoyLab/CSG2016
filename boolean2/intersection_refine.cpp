#include "precompile.h"

#include <map>
#include <set>

#include <xstruct.hpp>
#include <xlogger.h>

#include "intersection.h"
#include "hybrid_geometry.h"

namespace Boolean
{

    namespace
    {
        IndexPair makePbiIndex(const PbiRep* rep)
        {
            IndexPair res;
            MakeIndex(rep->ends, res);
            return res;
        }

        struct LinOrderItem : public LinOrderObjDefaultItem
        {
            VertexIndex vertex_idx;
        };
    }


    void EdgeInsctData::refine(Triangle* triangle, int which_edge)
    {
        if (isRefined()) return;

        assert(checkOrientation());

        PlaneLine plane_rep_edge(triangle->supportingPlane(), triangle->boundingPlane(which_edge));
        if (triangle->edge(which_edge).faceOrientation(triangle) > 0)
        {
            plane_rep_edge.inverse();
        }

        // construct sorting array
        std::vector<LinOrderItem> seqs(points.size());
        auto vertex_itr = points.begin();
        for (int i = 0; i < points.size(); i++)
        {
            seqs[i].vertex_idx = vertex_itr->vertex_idx;
            MyVertex& vertex = xvertex(vertex_itr->vertex_idx);
            if (vertex.isPlaneRep())
            {
                seqs[i].plane = vertex_itr->plane_rep;
            }
            else
            {
                seqs[i].point = vertex.vertex_rep();
                seqs[i].plane = XPlane();
            }

            vertex_itr++;
        }

        //std::map<VertexIndex, XPlane> vertex_id_plane_dict;
        //for (auto& set : inscts)
        //{
        //    for (auto &pbi : set.second)
        //    {
        //        vertex_id_plane_dict[pbiItr->ends[0]] = pbiItr->pends[0];
        //        vertex_id_plane_dict[pbiItr->ends[1]] = pbiItr->pends[1];
        //    }
        //}

        //// assign each vertex a plane and the corresponding id
        //for (int i = 0; i < points.size(); i++)
        //{
        //    auto res0 = vertex_id_plane_dict.find(seqs[i].id);
        //    if (res0 == vertex_id_plane_dict.end())
        //    {
        //        // some vertex cannot be assigned a valid plane
        //        // therefore we should find it manually from Plane Triples
        //        auto &vRef = xvertex(seqs[i].id);
        //        if (vRef.isPlaneRep())
        //        {
        //            auto& xpointRef = vRef.ppoint();
        //            for (int j = 0; j < 3; j++)
        //            {
        //                Real fres = data.line.dot(xpointRef.plane(j));
        //                if (fres == Real(0)) continue;

        //                if (fres > 0)
        //                    seqs[i].plane = xpointRef.plane(j);
        //                else if (fres < 0)
        //                    seqs[i].plane = xpointRef.plane(j).opposite();
        //                break;
        //            }
        //        }
        //    }
        //    else seqs[i].plane = res0->second;
        //}

        LinOrderObj<LinOrderItem> orderObj(plane_rep_edge);
        std::sort(seqs.begin(), seqs.end(), orderObj);

        // construct new pbis
        PbiList new_pbis(seqs.size()+1);
        std::map<VertexIndex, PbiList::size_type> vertex_vecidx_dict;
        const MyEdge& edge = triangle->edge(which_edge);
        vertex_vecidx_dict[edge.ends[0]] = 0;
        vertex_vecidx_dict[edge.ends[1]] = points.size() + 1;
        for (int i = 0; i < seqs.size(); i++)
        {
            vertex_vecidx_dict[seqs[i].vertex_idx] = i + 1;
            new_pbis[i].ends[1] = seqs[i].vertex_idx;
            new_pbis[i + 1].ends[0] = seqs[i].vertex_idx;
        }

        new_pbis[0].ends[0] = edge.ends[0];
        new_pbis.back().ends[1] = edge.ends[1];

        //for (int i = 0; i < seqs.size(); i++)
        //{
        //    vertex_vecidx_dict[seqs[i].id] = i + 1;
        //    newPbi[i].ends[1] = seqs[i].id;
        //    newPbi[i + 1].ends[0] = seqs[i].id;
        //}
        //newPbi[0].ends[0] = data.start;
        //newPbi[points.size()].ends[1] = data.end;

        // add pbi intersections' neighborInfo into newly constructed pbi
        for (auto pbi_itr = pbi_begin(); pbi_itr; ++pbi_itr)
        {
            assert(vertex_vecidx_dict.find(pbi_itr->ends[0]) != vertex_vecidx_dict.end());
            assert(vertex_vecidx_dict.find(pbi_itr->ends[1]) != vertex_vecidx_dict.end());

            int start = vertex_vecidx_dict[pbi_itr->ends[0]];
            int last = vertex_vecidx_dict[pbi_itr->ends[1]];

            assert(start < last);
            for (int i = start; i < last; i++)
            {
                assert(pbi_itr->neighbor.size() == 1);
                //new_pbis[i].neighbor.push_back(*pbi_itr->neighbor.begin());
                new_pbis[i].neighbor.insert(*pbi_itr->neighbor.begin());
            }
        }
        
        inscts.clear();
        inscts[INVALID_UINT32].swap(new_pbis);

        bRefined = true;
    }

    namespace
    {
        struct ResolveAuxInfo
        {
            XPlane        plane_a0, plane_a1, plane_b0, plane_b1;
            VertexIndex   a0, a1, b0, b1;
        };

        void fillAuxInfo(ResolveAuxInfo& info, FacePbi* pbi_a, FacePbi* pbi_b)
        {
            info.a0 = pbi_a->ends[0];
            info.a1 = pbi_a->ends[1];

            info.b0 = pbi_b->ends[0];
            info.b1 = pbi_b->ends[1];

            info.plane_a0 = pbi_a->pends[0];
            info.plane_a1 = pbi_a->pends[1];

            info.plane_b0 = pbi_b->pends[0];
            info.plane_b1 = pbi_b->pends[1];
        }

        class InsctResolveObj :
            public IResolveLineIntersection<
                MyVertex, MyVertex,
                MyVertex, MyVertex>
        {
        public:
            InsctResolveObj(std::array<LinOrderItem, 4>& item, std::array<int, 2>& counts,
                const XPlane& splane, FaceInsctData* insct, FacePbi* pbi_a, FacePbi* pbi_b) :
                item_(&item), counts_(&counts), splane_(splane), 
                insct_(insct), pbis_{ pbi_a, pbi_b } {}

            void resolve(
                Oriented_side side[4],
                const XPlane& plane_a, const XPlane& plane_b,
                const MyVertex& a0, const MyVertex& a1,
                const MyVertex& b0, const MyVertex& b1)
            {
                ResolveAuxInfo aux;
                fillAuxInfo(aux, pbis_[0], pbis_[1]);

                if (side[0] * side[1] == 0)
                {
                    if (side[2] * side[3] == 0)
                    {
                        // 相交于某个顶点，不用split
                        return;
                    }
                    assert(side[2] * side[3] == -1);

                    /// put point of a into b
                    LinOrderItem& item = item_->at(2 + counts_->at(1)++);
                    item.vertex_idx = (side[0] == ON_ORIENTED_BOUNDARY ? aux.a0 : aux.a1);
                    item.plane = plane_a;

                    PlaneLine line(splane_, plane_b);
                    make_good(line, xvertex(item.vertex_idx), &item);
                }
                else
                {
                    assert(side[0] * side[1] == -1);
                    if (side[2] * side[3] == 0)
                    {
                        LinOrderItem& item = item_->at(counts_->at(0)++);
                        item.vertex_idx = (side[2] == ON_ORIENTED_BOUNDARY ? aux.b0 : aux.b1);
                        item.plane = plane_b;

                        PlaneLine line(splane_, plane_a);
                        make_good(line, xvertex(item.vertex_idx), &item);
                    }
                    else
                    {
                        assert(side[2] * side[3] == -1);
                        // new vertex
                        XPlane thirdPlane = plane_b;
                        PlaneLine(splane_, plane_a).make_positive(thirdPlane); // 似乎不需要，可以尝试注释这一句
                        PlanePoint newPoint(splane_, plane_a, thirdPlane);

                        std::vector<VertexIndex*> vecs;
                        FaceInsctData::Vertex* face_new_pos = insct_->point(newPoint, TRIPLE_CROSS_0);
                        vecs.push_back(&face_new_pos->vId);

                        // 去所有的邻居看一看
                        VertexIndex minVal = face_new_pos->vId,
                            *cur_vertex_id = nullptr;

                        for (int i = 0; i < 2; i++)
                        {
                            int i2 = (i + 1) % 2;
                            assert(pbis_[i]->neighbor.size());
                            for (auto& neiInfo : pbis_[i]->neighbor)
                            {
                                if (neiInfo.second.type == NeighborInfo::Edge)
                                {
                                    MyEdge& cur_con = xedge(neiInfo.first);
                                    assert(cur_con.inscts);
                                    cur_vertex_id = cur_con.inscts->find_point(
                                        newPoint);
                                    
                                    // because it is edge-face intersection, 
                                    // it should have been detected in previous test
                                    if (!cur_vertex_id)
                                    {
                                        XLOG_ERROR << "There should have a intersection.";
                                        throw 1;
                                    }
                                }
                                else
                                {
                                    assert(neiInfo.second.type == NeighborInfo::Face);
                                    assert(neiInfo.second.pTrangle->inscts);

                                    face_new_pos = neiInfo.second.pTrangle->inscts->point(newPoint, TRIPLE_CROSS_0);
                                    cur_vertex_id = &face_new_pos->vId;
                                }

                                if (*cur_vertex_id < minVal)
                                {
                                    minVal = *cur_vertex_id;
                                }
                                vecs.push_back(cur_vertex_id);
                            }
                        }

                        VertexIndex this_vertex_idx = INVALID_UINT32;
                        if (minVal == INVALID_UINT32)
                        {
                            // 所有的邻居都没有，那就真没有了
                            this_vertex_idx = GlobalData::getObject()->insertVertex(newPoint);
                            for (VertexIndex *int_value : vecs)
                            {
                                *int_value = this_vertex_idx;
                            }
                        }
                        else
                        {
                            this_vertex_idx = minVal;
                            std::set<VertexIndex> coincident_vertex_ids;
                            for (VertexIndex *int_value : vecs)
                            {
                                if (*int_value == INVALID_UINT32)
                                {
                                    *int_value = minVal;
                                }
                                else
                                {
                                    coincident_vertex_ids.insert(*int_value);
                                }
                            }
                            mergeVertices(coincident_vertex_ids);
                            this_vertex_idx = *vecs[0];
                        }

                        PlaneLine line;

                        LinOrderItem& item = item_->at(counts_->at(0)++);
                        item.vertex_idx = this_vertex_idx;
                        item.plane = plane_b;
                        line = PlaneLine(splane_, plane_a);
                        make_good(line, xvertex(item.vertex_idx), &item);

                        LinOrderItem& item1 = item_->at(counts_->at(1)++);
                        item1.vertex_idx = this_vertex_idx;
                        item1.plane = plane_a;
                        line = PlaneLine(splane_, plane_b);
                        make_good(line, xvertex(item1.vertex_idx), &item1);
                    }
                }
            }

            void make_good(const PlaneLine& line, const MyVertex& vertex, LinOrderItem* item)
            {
                if (!vertex.isPlaneRep())
                {
                    item->plane = XPlane();
                    item->point = vertex.vertex_rep();
                }
                else
                {
                    line.make_positive(item->plane);
                }
            }

        private:
            std::array<LinOrderItem, 4>*    item_;
            std::array<int, 2>*             counts_;
            XPlane                          splane_;
            FaceInsctData*                  insct_;
            FacePbi*                        pbis_[2];
        };

        class CoplanarResolveObj : public IResolveLineIntersection<
            MyVertex, MyVertex, MyVertex, MyVertex>
        {
        public:
            CoplanarResolveObj(std::array<LinOrderItem, 4>& item, std::array<int, 2>& counts,
                const XPlane& splane, FacePbi* pbi_a, FacePbi* pbi_b) :
                item_(&item), counts_(&counts), splane_(splane), pbis_{pbi_a, pbi_b} {}

            void resolve(Oriented_side side[4],
                const XPlane& plane_a, const XPlane& plane_b,
                const MyVertex& a0, const MyVertex& a1,
                const MyVertex& b0, const MyVertex& b1)
            {
                ResolveAuxInfo aux;
                fillAuxInfo(aux, pbis_[0], pbis_[1]);

                PlaneLine line(splane_, plane_a);
                Real dotRes = line.dot(aux.plane_b0);
                assert(dotRes != 0);
                bool inverseLine = dotRes > 0 ? false : true;

                if (inverseLine)
                {
                    aux.plane_b0.inverse();
                    aux.plane_b1.inverse();
                    std::swap(aux.plane_b0, aux.plane_b1);
                }

                assert(check_orientation(line, aux));

                // linear order 计算overlap
                if (line.linear_order_unsafe(aux.plane_b0, aux.plane_a1) <= 0 ||
                    line.linear_order_unsafe(aux.plane_a0, aux.plane_b1) <= 0)
                {
                    return;
                }

                int compRes[2] = {
                    line.linear_order_unsafe(aux.plane_a0, aux.plane_b0),
                    line.linear_order_unsafe(aux.plane_a1, aux.plane_b1),
                };

                // must not be coincident
                assert(!(compRes[0] == 0 && compRes[1] == 0));

                // ------ a
                //   ---- b
                if (compRes[0] > 0)
                {
                    LinOrderItem& item = item_->at(counts_->at(0)++);
                    item.vertex_idx = aux.b0;
                    if (b0.isPlaneRep())
                    {
                        item.plane = aux.plane_b0;
                    }
                    else
                    {
                        item.plane = XPlane();
                        item.point = b0.vertex_rep();
                    }
                }
                else if (compRes[0] < 0)
                {
                    LinOrderItem& item = item_->at(2+counts_->at(1)++);
                    item.vertex_idx = aux.a0;
                    if (a0.isPlaneRep())
                    {
                        item.plane = aux.plane_a0;
                        item.plane.inverse();
                    }
                    else
                    {
                        item.plane = XPlane();
                        item.point = a0.vertex_rep();
                    }
                }

                // ----   a
                // ------ b
                if (compRes[1] > 0)
                {
                    LinOrderItem& item = item_->at(2+counts_->at(1)++);
                    item.vertex_idx = aux.a1;
                    if (a1.isPlaneRep())
                    {
                        item.plane = aux.plane_a1;
                        item.plane.inverse();
                    }
                    else
                    {
                        item.plane = XPlane();
                        item.point = a1.vertex_rep();
                    }
                }
                else if (compRes[1] < 0)
                {
                    LinOrderItem& item = item_->at(counts_->at(0)++);
                    item.vertex_idx = aux.b1;
                    if (b1.isPlaneRep())
                    {
                        item.plane = aux.plane_b1;
                    }
                    else
                    {
                        item.plane = XPlane();
                        item.point = b1.vertex_rep();
                    }
                }

            }

        private:
            bool check_orientation(const PlaneLine& line, ResolveAuxInfo& aux) const
            {
                if ((line.dot(aux.plane_a0) > 0) &&
                    (line.dot(aux.plane_a1) > 0) &&
                    (line.dot(aux.plane_b0) > 0) &&
                    (line.dot(aux.plane_b1) > 0) &&
                    (line.linear_order_unsafe(aux.plane_a0, aux.plane_a1) > 0) &&
                    (line.linear_order_unsafe(aux.plane_b0, aux.plane_b1) > 0))
                    return true;
                return false;
            }

            std::array<LinOrderItem, 4>*    item_;
            std::array<int, 2>*             counts_;
            XPlane                          splane_;
            FacePbi*                        pbis_[2];
        };
    }

    void FaceInsctData::resolveIntersection(Triangle* triangle)
    //void FaceInsctData::resolveIntersection(Triangle* pTri, std::vector<VertexIndex>* strayVertices)
    {
        if (inscts.size() <= 1) return;

        struct FacePbiTessData
        {
            PbiList::iterator pbi_itr;
            PbiList* container_itr;
            std::vector<LinOrderItem> points;
            //MeshIndex mesh_idx;
        };

        typedef PbiLists::iterator PbiListsItr;
        typedef PbiList::iterator PbiListItr;

        std::map<IndexPair, FacePbiTessData*> tessData;
        for (PbiListsItr setItr = inscts.begin(); setItr != inscts.end(); ++setItr)
        {
            PbiListsItr setItr2 = setItr; ++setItr2;
            for (; setItr2 != inscts.end(); ++setItr2)
            {
                for (PbiListItr pbiItr = setItr->second.begin(); pbiItr != setItr->second.end(); ++pbiItr)
                {
                    for (PbiListItr pbiItr2 = setItr2->second.begin(); pbiItr2 != setItr2->second.end(); ++pbiItr2)
                    {
                        PbiListItr pbi_itrs[2] = { pbiItr, pbiItr2 };
                        PbiListsItr set_itrs[2] = { setItr, setItr2 };

                        std::array<LinOrderItem, 4> insct_item_slots;
                        std::array<int,2> insct_count{ 0, 0 };

                        InsctResolveObj insct_obj(insct_item_slots, insct_count, 
                            triangle->supportingPlane(), this, &*pbiItr, &*pbiItr2);

                        CoplanarResolveObj coplanar_obj(insct_item_slots, insct_count,
                            triangle->supportingPlane(), &*pbiItr, &*pbiItr2);

                        LineInsctResultType insct_res = plane_based_line_intersection(
                            pbiItr->vertPlane, pbiItr2->vertPlane,
                            xvertex(pbi_itrs[0]->ends[0]), xvertex(pbi_itrs[0]->ends[1]),
                            xvertex(pbi_itrs[1]->ends[0]), xvertex(pbi_itrs[1]->ends[1]),
                            &insct_obj, &coplanar_obj
                        );

                        if (insct_res == No_Intersect) continue;
                        assert(insct_count[0] <= 2 && insct_count[1] <= 2);

                        for (int i = 0; i < 2; ++i)
                        {
                            if (insct_count[i] == 0) continue;

                            IndexPair pbiIndex = makePbiIndex(&*pbi_itrs[i]);
                            auto pbi_search_result = tessData.find(pbiIndex);
                            if (pbi_search_result == tessData.end())
                            {
                                FacePbiTessData *tessItem = new FacePbiTessData;
                                tessItem->container_itr = &(set_itrs[i]->second);
                                tessItem->pbi_itr = pbi_itrs[i];
                                //tessItem->mesh_idx = set_itrs[i]->first;

                                auto insertRes = tessData.insert(
                                    decltype(tessData)::value_type(
                                        pbiIndex, tessItem
                                    )
                                );

                                assert(insertRes.second); // insertion succeed
                                pbi_search_result = insertRes.first;
                            }

                            for (int j = 0; j < insct_count[i]; ++j)
                            {
                                pbi_search_result->second->points.
                                    push_back(insct_item_slots[2*i+j]);
                            }
                        }

                    } // loop 4
                } // loop 3
            } // loop 2
        } // outer loop

        std::map<PbiList*, std::vector<int>> garbage;
        for (auto &item : tessData)
        {
            FacePbiTessData* tess_pbi_info = item.second;
            PlaneLine line(triangle->supportingPlane(), tess_pbi_info->pbi_itr->vertPlane);

            LinOrderObj<LinOrderItem> orderObj(line);
            auto &inserted = tess_pbi_info->points;
            std::sort(inserted.begin(), inserted.end(), orderObj);

            // coincident vertex can be generated from '-|-' situation
            inserted.erase(
                std::unique(
                    inserted.begin(),
                    inserted.end(),
                    [](const LinOrderItem &a, const LinOrderItem &b)->bool {
                        return a.vertex_idx == b.vertex_idx;
                    }
                ),
                inserted.end()
            );

            // template inherent prep, neighbor
            FacePbi template_insct_data = *tess_pbi_info->pbi_itr;
            template_insct_data.pends[0] = XPlane();
            template_insct_data.pends[1] = XPlane();

            std::vector<FacePbi> newPbi(inserted.size() + 1, template_insct_data);
            std::map<VertexIndex, uint32_t> vertex_vecidx_dict;
            vertex_vecidx_dict[tess_pbi_info->pbi_itr->ends[0]] = 0;
            vertex_vecidx_dict[tess_pbi_info->pbi_itr->ends[1]] = inserted.size() + 1;
            for (int i = 0; i < inserted.size(); i++)
            {
                vertex_vecidx_dict[inserted[i].vertex_idx] = i + 1;
                newPbi[i].ends[1] = inserted[i].vertex_idx;
                newPbi[i + 1].ends[0] = inserted[i].vertex_idx;
            }
            newPbi[0].ends[0] = tess_pbi_info->pbi_itr->ends[0];
            newPbi.back().ends[1] = tess_pbi_info->pbi_itr->ends[1];

            PbiList* pbi_list = tess_pbi_info->container_itr;
            //XR::vec_quick_delete(tess_pbi_info->pbi_itr, *pbi_list);
            garbage[pbi_list].push_back(tess_pbi_info->pbi_itr-pbi_list->begin());
            pbi_list->insert(pbi_list->end(), newPbi.begin(), newPbi.end());
            newPbi.clear();
        }

        for (auto& group : garbage)
        {
            std::sort(group.second.begin(), group.second.end(), std::greater<int>());
            for (auto pbi_itr : group.second)
            {
                XR::vec_quick_delete(pbi_itr, *group.first);
            }
        }
    }

    void FaceInsctData::removeOverlapPbi()
    {
        typedef PbiList::iterator PbiIterator;
        std::map<VertexIndex, std::set<FacePbi*>> sparse_graph;
        //std::vector<decltype(inscts[0].begin())> garbage;

        for (auto &insct_itr : inscts)
        {
            auto& pbi_list = insct_itr.second;
            std::vector<PbiIterator> garbage;
            for (PbiIterator pbi_itr = pbi_list.begin(); pbi_itr != pbi_list.end(); ++pbi_itr)
            {
                bool found = false;
                for (FacePbi* alreadyHere : sparse_graph[pbi_itr->ends[0]]) // search in current graph
                {
                    if (alreadyHere->ends[0] == pbi_itr->ends[1] ||
                        alreadyHere->ends[1] == pbi_itr->ends[1]) // if has
                    {
                        for (auto& nInfo : pbi_itr->neighbor)
                        {
                            bool unique = true; // search if the neighInfo is already there, == std::find
                            for (auto& nInfo2 : alreadyHere->neighbor)
                            {
                                if (nInfo.first == nInfo.first)
                                {
                                    unique = false;
                                    break;
                                }
                            }

                            if (unique)
                            {
                                //alreadyHere->neighbor.push_back(nInfo);
                                alreadyHere->neighbor.insert(nInfo);
                            }
                        }
                        found = true;
                        break;
                    }
                }

                if (found)
                {
                    garbage.push_back(pbi_itr);
                    continue;
                }
                else
                {
                    sparse_graph[pbi_itr->ends[0]].insert(&*pbi_itr);
                    sparse_graph[pbi_itr->ends[1]].insert(&*pbi_itr);
                }
            }

            const int number_of_overlap = garbage.size();
            std::sort(garbage.begin(), garbage.end(), std::greater<PbiIterator>());
            for (PbiIterator pbi_itr : garbage)
            {
                XR::vec_quick_delete(pbi_itr, pbi_list);
            }
        }


    }

    void FaceInsctData::refine(Triangle* triangle)
    {
        if (isRefined()) return;

        assert(checkOrientation());

        removeOverlapPbi();
        if (inscts.size() >= 2)
        {
            resolveIntersection(triangle);
            removeOverlapPbi();
        }
        bRefined = true;
    }
}