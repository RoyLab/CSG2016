#include "precompile.h"
#include <vector>
#include <memory>
#include <queue>
#include "RegularMesh.h"
#include "Octree.h"
#include "csg.h"

namespace Boolean
{
    template <class T> struct AutoPtr : std::shared_ptr<T> {};
    const int MARK_BEGIN = 0xff; // 因为mark还用来在第三阶段标志有没有被访问过，所以这里让出256个数字用于这些工作
    enum Mark { UNVISITED, SEEDED0, SEEDED1, SEEDED2, VISITED };

    struct SSeed
    {
        MyEdge::Index edgeId;
        IPolygon* pFace = nullptr;
        AutoPtr<IIndicatorVector> eIndicators;

        SSeed() {}
        SSeed(const SSeed& other) { *this = other; }
        SSeed& operator=(const SSeed& other);
    };

    void doClassification(Octree* pOctree, CSGTree<RegularMesh>* pCSG,
        std::vector<RegularMesh*>& meshList, RegularMesh* result,
        MyVertex::Index seedId)
    {
        typedef uint32_t MeshId;

        uint32_t nMesh = meshList.size();
        CSGTreeOld* tree = pCSG->auxiliary();
        CSGTreeNode** curTreeLeaves = new CSGTreeNode*[nMesh];

        SSeed tmpSeed;
        MyVertex& seedV = xvertex(seedId);
        tmpSeed.edgeId = *seedV.edges.begin();
        MyEdge& seedE = xedge(tmpSeed.edgeId);
        tmpSeed.pFace = MyEdge::FaceIterator(seedE).face();

        tmpSeed.eIndicators.reset(new FullIndicatorVector(nMesh));
        calcEdgeIndicatorByExtremity(seedV, tmpSeed, tmpSeed.eIndicators.get());

        MeshId curMeshId;
        std::queue<IPolygon*> faceQueue;
        std::queue<SSeed> intraQueue, interQueue;
        bool added, inverse;
        Relation relation;
        TestTree dummyForest;
        std::vector<Relation> relTab(nMesh);
        std::vector<MyEdge::Index> edges;

        interQueue.push(tmpSeed);
        tmpSeed.pFace->mark = SEEDED0;
        while (!interQueue.empty())
        {
            SSeed curSeed = interQueue.front();
            interQueue.pop();

            if (curSeed.pFace->mark == VISITED)
                continue;

            curMeshId = tmpSeed.pFace->meshId();
            assert(intraQueue.empty());
            curSeed.pFace->mark = SEEDED1;
            intraQueue.push(curSeed);

            while (!intraQueue.empty())
            {
                curSeed = intraQueue.front();
                intraQueue.pop();

                if (curSeed.pFace->mark == VISITED)
                    continue;

                CSGTreeNode* tree0 = copy2(tree->pRoot, curTreeLeaves);
                calcFaceIndicator(curSeed, relTab);

                relation = ParsingCSGTree(meshList[curMeshId], relTab.data(), 
                    nMesh, tree0, curTreeLeaves, dummyForest);
                assert(relation != REL_NOT_AVAILABLE || relation != REL_UNKNOWN);

                inverse = meshList[curMeshId]->inverse();
                added = (relation == REL_SAME);

                faceQueue.push(curSeed.pFace);
                while (!faceQueue.empty())
                {
                    IPolygon* curFace = faceQueue.front();
                    faceQueue.pop();
                    if (curFace->mark == VISITED)
                        continue;

                    curFace->mark = VISITED;
                    if (added)
                    {
                        result->faces().push_back(curFace);
                        result->inverseMap.push_back(inverse);
                    }

                    edges.clear();
                    curFace->getEdges(edges);
                    assert(edges.size() == curFace->degree());
                    for (int i = 0; i < curFace->degree(); i++)
                    {
                        MyEdge& curEdge = xedge(edges[i]);
                        MyEdge::FaceIterator fItr(curEdge);
                        if (curEdge.neighbor)
                        {
                            assert(curEdge.faceCount() >= 4);
                            for (; fItr; ++fItr)
                            {
                                tmpSeed.edgeId = edges[i];
                                tmpSeed.pFace = fItr.face();
                                if (fItr.face()->meshId() == curMeshId)
                                {
                                    if (fItr.face()->mark < SEEDED1)
                                    {
                                        tmpSeed.pFace->mark = SEEDED1;
                                        intraQueue.push(tmpSeed);
                                    }
                                }
                                else
                                {
                                    if (fItr.face()->mark < SEEDED0)
                                    {
                                        tmpSeed.pFace->mark = SEEDED0;
                                        interQueue.push(tmpSeed);
                                    }
                                }
                            }
                        }
                        else
                        {
                            for (; fItr; ++fItr)
                            {
                                if (fItr.face()->mark < SEEDED2)
                                {
                                    faceQueue.push(fItr.face());
                                    fItr.face()->mark == SEEDED2;
                                }
                            }
                        }
                    }
                }
            }
        }
        SAFE_DELETE_ARRAY(curTreeLeaves);
    }

    SSeed & SSeed::operator=(const SSeed & other)
    {
        edgeId = other.edgeId;
        pFace = other.pFace;

        if (other.eIndicators->getType() == IIndicatorVector::FULL)
            eIndicators.reset(new FullIndicatorVector(
                *reinterpret_cast<FullIndicatorVector*>(other.eIndicators.get())));
        else eIndicators.reset(new SampleIndicatorVector(
                *reinterpret_cast<SampleIndicatorVector*>(other.eIndicators.get())));
    }

}