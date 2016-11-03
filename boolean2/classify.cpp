#include "precompile.h"
#include <vector>
#include "RegularMesh.h"
#include "Octree.h"

namespace Boolean
{
    void doClassification(Octree* pOctree, CSGTree<RegularMesh>* pCSG, std::vector<RegularMesh*>& meshList, RegularMesh* result, uint32_t* extremes)
    {
        assert(extremes);
        MyVertex::Index seedId = extremes[0];

        MyVertex& seedV = xvertex(seedId);
        MyEdge& seedE = xedge(*seedV.edges.begin());

        size_t nMesh = meshList.size();
        pConstruct = new Delegate<HalfedgeDS>;

        GroupParseInfo infos;

        auto tree = pCsg->auxiliary();
        infos.curTreeLeaves = new CSGTreeNode*[nMesh];

        for (size_t imesh = 0; imesh < nMesh; imesh++)
        {
            SeedInfoWithId idSeed;
            idSeed.meshId = imesh;
            createFirstSeed(idSeed);
            infos.otherMeshSeedQueue.push(idSeed);

            while (!infos.otherMeshSeedQueue.empty())
            {
                SeedInfoWithId& curSeed = infos.otherMeshSeedQueue.front();
                std::vector<int> itstPrims;
                itstAlg->get_adjGraph()->getIntersectPrimitives(curSeed.meshId, itstPrims);

                CSGTreeNode* tree0 = copy2(tree->pRoot, infos.curTreeLeaves);
                for (size_t id = 0; id < nMesh; id++)
                    infos.curTreeLeaves[id]->relation = static_cast<Relation>(curSeed.indicators->at(id));

                for (size_t id : itstPrims)
                    infos.curTreeLeaves[id]->relation = REL_UNKNOWN;

                Relation meshRel = CompressCSGNodeIteration(tree0);
                infos.ttree1.reset(tree0);

                if (meshRel == REL_NOT_AVAILABLE)
                {
                    infos.curMeshId = curSeed.meshId;
                    infos.curMeshSeedQueue.push(curSeed);
                    curSeed.seedFacet->mark = SEEDED1;

                    while (!infos.curMeshSeedQueue.empty())
                    {
                        SeedInfo& sndInfo = infos.curMeshSeedQueue.front();
                        if (sndInfo.seedFacet->mark != VISITED)
                        {
                            if (sndInfo.seedFacet->isSimple())
                                floodSimpleGroup(infos, sndInfo);
                            else
                                floodComplexGroup(infos, sndInfo);
                        }
                        infos.curMeshSeedQueue.pop();
                    }
                }
                infos.otherMeshSeedQueue.pop();
            }
        }
        csgResult->delegate(*pConstruct);
        SAFE_DELETE(pConstruct);
        SAFE_DELETE_ARRAY(infos.curTreeLeaves);


    }

}