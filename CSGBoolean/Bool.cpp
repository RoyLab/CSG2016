#define CSG_EXPORTS
#include "bool.h"
#include "MyAlgorithm.h"
#include "MyMesh.h"
#include <string>
//#include "precompile.h"
//#include <OpenMesh/Core/IO/MeshIO.hh>
//#include "MyMesh.h"
//#include <ctime>
//#include <list>
//#include "COctree.h"
//#include "BinaryTree.h"
//#include "isect.h"
//#include "IsectTriangle.h"
//#include <queue>
//#include "CSGExprNode.h"
//#include "configure.h"
//#include <random>
//#include <ctime>
//#include <OpenMesh/Tools/Utils/getopt.h>

//#ifdef _DEBUG
////#include <vld.h>
//int countd1, countd2, countd3, countd4, countd5;
//#endif


//#define DEBUG_ISECT
//#define DEBUG_FLOOD
//#define FIXED_COLOR

//#pragma warning(disable: 4800 4996)


namespace CSG
{
    extern "C" CSG_API int loadMesh(std::vector<MyMesh*>& meshes, const std::vector<std::string>& names)
    {
        return 0;
    }


    void test()
    {
        std::vector<std::string> names;
        std::string expr;

        std::vector<MyMesh*> meshList;
        loadMesh(meshList, names);

        for (auto pMesh : meshList)
            pMesh->update();

        MyAlgorithm *alg = new MyAlgorithm;
        alg->solve(expr, meshList);

        MyMesh* res = alg->popResultMesh();
        SAFE_DELETE(alg);
    }


//    extern ISectZone* ZONE;
//
//    HANDLE _output;
//    clock_t t0;
//
//    void GetLeafNodes(OctreeNode* pNode, std::list<OctreeNode*>& leaves, int NodeType);
//
//    void StdOutput(const char* str)
//    {
//        std::string ch(str);
//        WriteConsole(_output, str, ch.size(), 0, 0);
//        WriteConsole(_output, "\n", 1, 0, 0);
//    }
//
//    void DebugInfo(char* str, clock_t& t0)
//    {
//        char ch[32];
//        sprintf(ch, "%s: %d\0", str, clock()-t0);
//        t0 = clock();
//        StdOutput(ch);
//    }
//
//
//    static void RelationTest(OctreeNode* pNode, Octree<>* pOctree, std::map<unsigned, Relation>& famap)
//    {
//        for (auto &pair: pNode->DiffMeshIndex) pair.Rela = famap[pair.ID];
//
//        std::map<unsigned, Relation> relmap;
//        if (pNode->Child)
//        {
//            for (int i = 0; i < 8; i++)
//            {
//                for (auto &pair: pNode->Child[i].DiffMeshIndex)
//                {
//                    if (relmap.find(pair.ID) != relmap.end()) continue;
//                    else relmap[pair.ID] = PolyhedralInclusionTest(pNode->Child[i].BoundingBox.Center(), pOctree, pair.ID);
//                }
//            }
//
//            for (unsigned i = 0; i < 8; i++)
//                RelationTest(&pNode->Child[i], pOctree, relmap);
//        }
//    }
//
//    static void RelationTest(Octree<>* pOctree)
//    {
//        std::map<unsigned, Relation> relmap;
//        RelationTest(pOctree->Root, pOctree, relmap);
//    }
//
//    void ISectTest(Octree<>* pOctree)
//    {
//        assert(pOctree);
//        std::list<OctreeNode*> leaves;
//        GetLeafNodes(pOctree->Root, leaves, NODE_COMPOUND);
//
//        std::map<GS::IndexPair, std::set<GS::IndexPair>> antiOverlapMap;
//        for (auto leaf: leaves)
//        {
//            auto itr = leaf->TriangleTable.begin();
//            auto iEnd = leaf->TriangleTable.end();
//            decltype(iEnd) itr2;
//            unsigned i, j, ni, nj;
//            MyMesh *meshi, *meshj;
//            MyMesh::FaceHandle tri1, tri2;
//            Vec3d *v0,*v1,*v2, nv,*u0,*u1,*u2,nu,start,end;
//            MyMesh::FVIter fvItr;
//            int isISect;
//            int startT(0), endT(0);
//            ISectTriangle **si=nullptr, **sj=nullptr;
//            VertexPos startiT, startjT, endiT, endjT;
//            ISVertexItr vP1, vP2;
//
//            int meshId[2];
//            for (; itr != iEnd; ++itr)
//            {
//                itr2 = itr;
//                ++itr2;
//                for (; itr2 != iEnd; ++itr2)
//                {
//                    ni = itr->second.size();
//                    nj = itr2->second.size();
//
//                    meshi = pOctree->pMesh[itr->first];
//                    meshj = pOctree->pMesh[itr2->first];
//
//                    //if (meshi->ID != 100 && meshj->ID != 100) continue;
//
//                    if (meshi->ID > meshj->ID) {meshId[0] = meshj->ID; meshId[1] = meshi->ID;}
//                    else {meshId[0] = meshi->ID; meshId[1] = meshj->ID;}
//
//                    GS::IndexPair iPair;
//                    GS::MakeIndex(meshId, iPair);
//                    auto &antiOverlapSet = antiOverlapMap[iPair];
//
//                    for (i = 0; i < ni; i++)
//                    {
//                        for (j = 0; j < nj; j++)
//                        {
//                            tri1 = itr->second[i];
//                            tri2 = itr2->second[j];
//                            if (meshi->ID > meshj->ID) {meshId[0] = tri2.idx(); meshId[1] = tri1.idx();}
//                            else {meshId[0] = tri1.idx(); meshId[1] = tri2.idx();}
//                            GS::MakeIndex(meshId, iPair);
//                            if (antiOverlapSet.find(iPair) != antiOverlapSet.end()) continue;
//                            else antiOverlapSet.insert(iPair);
//
//                            // intersection test main body
//                            GetCorners(meshi, tri1, v0, v1, v2);
//                            GetCorners(meshj, tri2, u0, u1, u2);
//
//                            nv = meshi->normal(tri1);
//                            nu = meshj->normal(tri2);
//                            
//                            startT = INNER; endT = INNER; // return to Zero.
//
//                            //auto &isecTris = (*si)->isecTris[meshj->ID];
//                            //auto &icoplTris = (*si)->coplanarTris;
//
//                            isISect = TriTriIntersectTest(*v0, *v1, *v2, nv,
//                                *u0, *u1, *u2, nu, startT, endT, start, end);
//
//                            if (isISect < 0) continue;
//
//                            si = &meshi->property(meshi->SurfacePropHandle, tri1);
//                            sj = &meshj->property(meshj->SurfacePropHandle, tri2);
//
//                            if (!*si) *si = new ISectTriangle(meshi, tri1);
//                            if (!*sj) *sj = new ISectTriangle(meshj, tri2);
//
//                            if (isISect == 0)
//                            {
//                                (*si)->coplanarTris[meshj->ID].emplace_back(tri2);
//                                (*sj)->coplanarTris[meshi->ID].emplace_back(tri1);
//                                continue;
//                            }
//
//                            startiT = VertexPos(startT & 0xffff);
//                            startjT = VertexPos(startT >> 16);
//
//                            endiT = VertexPos(endT & 0xffff);
//                            endjT = VertexPos(endT >> 16);
//
//                            if (IsEqual(start, end))
//                            {
//                                // 点相交
//                                Vec3d point = (start+end)/2;
//
//                                // 最后一个参数表示，可能存在两个以上的插入点
//                                vP1 = InsertPoint(*si, startiT, point);
//                                InsertPoint(*si, endiT, vP1);
//                                InsertPoint(*sj, startjT, vP1);
//                                InsertPoint(*sj, endjT, vP1);
//                            }
//                            else
//                            {
//                                // 线相交
//                                double d = OpenMesh::dot(OpenMesh::cross(nv, nu), end-start);
//                                if (meshi->bInverse ^ meshj->bInverse) d = -d;
//                                if (d > 0)
//                                {
//                                    vP1 = InsertPoint(*si, startiT, start);
//                                    vP2 = InsertPoint(*si, endiT, end);
//                                    InsertSegment(*si, vP1, vP2, *sj);
//                                    vP1 = InsertPoint(*sj, startjT, vP1);
//                                    vP2 = InsertPoint(*sj, endjT, vP2);
//                                    InsertSegment(*sj, vP2, vP1, *si);
//                                }
//                                else
//                                {
//                                    vP1 = InsertPoint(*si, startiT, start);
//                                    vP2 = InsertPoint(*si, endiT, end);
//                                    InsertSegment(*si, vP2, vP1, *sj);
//                                    vP1 = InsertPoint(*sj, startjT, vP1);
//                                    vP2 = InsertPoint(*sj, endjT, vP2);
//                                    InsertSegment(*sj, vP1, vP2, *si);
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//
//    void FloodColoring(Octree<>* pOctree, CSGTree* pPosCSG);
//
//    GS::BaseMesh* result;
//    int fDebugMode = -1;
//    unsigned startId = 0;
//    int percent = 100;
//
//    static GS::BaseMesh* BooleanOperation2(GS::CSGExprNode* input, HANDLE stdoutput, bool bd)
//    {
//        _output= stdoutput;
//        MyMesh** arrMesh = NULL;
//        result = new GS::BaseMesh;
//        int nMesh = -1;
//        StdOutput("Start:");
//        t0 = clock();
//        auto t = t0;
//        CSGTree* pCSGTree = ConvertCSGTree(input, &arrMesh, &nMesh);
//        CSGTree* pPosCSG = ConvertToPositiveTree(pCSGTree);
//        delete pCSGTree;
//        DebugInfo("Convert", t0);
//        Octree<>* pOctree = BuildOctree(arrMesh, nMesh);
//        DebugInfo("BuildTree", t0);
//        InitZone();
//        ISectTest(pOctree);
//        DebugInfo("ISectTest", t0);
//        if (bd)
//        {
//            std::ifstream dbFile("D:\\bool\\DebugConfig.ini");
//            if (dbFile)
//            {
//                dbFile >> fDebugMode;
//#ifdef DEBUG_FLOOD
//                dbFile >> percent;
//                dbFile >> startId;
//#endif
//            }
//            dbFile.close();
//        }
//        else fDebugMode = -1;
//        FloodColoring(pOctree, pPosCSG);
//        DebugInfo("FloodColoring", t0);
//
//        delete pOctree;
//        delete pPosCSG;
//        delete input;
//
//        ReleaseZone();
//        for (int i = 0; i < nMesh; i++)
//            delete arrMesh[i];
//        delete [] arrMesh;
//        auto k = clock()-t;
//        char ch[32];
//        wsprintf(ch, "time: %d\n", k);
//        StdOutput(ch);
//        return result;
//    }
//
//
//    extern "C" CSG_API GS::BaseMesh* BooleanOperation(GS::CSGExprNode* input, HANDLE stdoutput, bool db)
//    {
//        return BooleanOperation2(input, stdoutput, db);
//    }
//
//    extern "C" CSG_API GS::BaseMesh* BooleanOperation_MultiThread(GS::CSGExprNode* input)
//    {
//        return NULL;
//    }
//
//    extern "C" CSG_API void SnapModel(const GS::BaseMesh* mesh)
//    {
//        MyMesh *kk = ConvertToMPMeshChrome(mesh);
//        OpenMesh::IO::Options wopt;
//        //wopt += OpenMesh::IO::Options::FaceNormal;
//        wopt += OpenMesh::IO::Options::FaceColor;
//        if (OpenMesh::IO::write_mesh(*kk, "C:\\Users\\RUI\\Desktop\\model.obj", wopt))
//            StdOutput("Write files success.\n");
//        else StdOutput("Write files failed.\n");
//        kk->release_face_colors();
//        delete kk;
//    }
//
//
//    void GetLeafNodes(OctreeNode* pNode, std::list<OctreeNode*>& leaves, int NodeType)
//    {
//        if (pNode == NULL) return;
//
//        if (IsLeaf(pNode) || pNode->Type == NODE_SIMPLE)
//        {
//            if (!pNode->TriangleCount) return; 
//            if (NodeType == pNode->Type)
//                leaves.push_back(pNode);
//            return;
//        }
//
//        for (int i = 0; i < 8 ; i++)
//            GetLeafNodes(&pNode->Child[i], leaves, NodeType);
//    }
//
//    struct FacePair
//    {
//        MyMesh::FaceHandle seed, current;
//        const MyMesh::FaceHandle& operator[](int index) const
//        {
//            switch (index)
//            {
//            case 0:    return seed;
//            case 1:    return current;
//            default:
//                assert(0);
//                return current;
//            }
//        }
//
//        MyMesh::FaceHandle& operator[](int index)
//        {
//            switch (index)
//            {
//            case 0:    return seed;
//            case 1:    return current;
//            default:
//                assert(0);
//                return current;
//            }
//        }
//    };
//
//    struct SeedInfo
//    {
//        std::queue<FacePair> queue;
//        Relation *relation;
//
//        SeedInfo():relation(nullptr){}
//        ~SeedInfo(){SAFE_RELEASE_ARRAY(relation);}
//    };
//
//    bool CompareRelationSpace(ISectTriangle* triSeed, ISectTriangle* triCur)
//    {
//        assert(triSeed && triCur);
//
//        std::map<int, int> map;
//        for (auto& pair: triSeed->segs)
//            map[pair.first] = 1;
//        for (auto& pair2: triSeed->coplanarTris)
//            map[pair2.first] = 1;
//
//        std::map<int, int>::iterator itr;
//        for (auto& pair: triCur->segs)
//        {
//            itr = map.find(pair.first);
//            if (itr == map.end())
//                return false;
//            else itr->second ++;
//        }
//
//        for (auto& pair2: triCur->coplanarTris)
//        {
//            itr = map.find(pair2.first);
//            if (itr == map.end())
//                return false;
//            else itr->second ++;
//        }
//
//        for (auto& pair3: map)
//            if (pair3.second == 1) return false;
//
//        return true;
//    }
//
//
//        static void AddTriangle(MyMesh* pMesh, MyMesh::FaceHandle face)
//    {
//#ifdef _DEBUG
//        countd1 ++;
//#endif
//        GS::double3 v[3];
//        Vec3d *v0, *v1, *v2;
//        if (pMesh->bInverse) GetCorners(pMesh, face, v2, v1, v0);
//        else GetCorners(pMesh, face, v0, v1, v2);
//        v[0] = Vec3dToDouble3(*v0);
//        v[1] = Vec3dToDouble3(*v1);
//        v[2] = Vec3dToDouble3(*v2);
//        result->AddTriangle(v);
//    }
//
//    static void AddTriangle(MyMesh* pMesh, MyMesh::FaceHandle face, GS::float4& color)
//    {
//#ifdef _DEBUG
//        countd1 ++;
//#endif
//        GS::double3 v[3];
//        Vec3d *v0, *v1, *v2;
//        if (pMesh->bInverse) GetCorners(pMesh, face, v2, v1, v0);
//        else GetCorners(pMesh, face, v0, v1, v2);
//        v[0] = Vec3dToDouble3(*v0);
//        v[1] = Vec3dToDouble3(*v1);
//        v[2] = Vec3dToDouble3(*v2);
//        result->AddTriangle(v, color);
//    }
//
//
//    int offset = -5;
//
//
//    void FloodColoring(Octree<>* pOctree, CSGTree* pPosCSG)
//    {
//#ifdef FIXED_COLOR
//        std::ifstream cfile("D:\\bool\\color.ini");
//        if (!cfile)
//        {
//            StdOutput("Oh my God,0.0");
//            return;
//        }
//
//        GS::float4 colorBuffer[100];
//        int i = 100;
//        float r,g,b;
//        unsigned cstart = 0;
//        cfile >> cstart;
//        while (i--)
//        {
//            GS::float4& cc = colorBuffer[99-i];
//            cfile >> r >> g >> b;
//            cc.x = r; cc.y = g, cc.z = b;
//            cc.w = 1.0f;
//        }
//        unsigned icolor = cstart;
//        cfile.close();
//#endif
//        const int seeded = 4;
//
//        MyMesh *pMesh;
//        std::queue<MyMesh::FaceHandle> faceQueue;
//        std::queue<SeedInfo> seedQueueList;
//        MyMesh::FaceHandle curFace, relatedFace;
//        MyMesh::FaceFaceIter ffItr;
//        Vec3d *v0, *v1, *v2;
//        Relation *curRelationTable, curRelation;
//        CSGTreeNode* curTree;
//        CSGTreeNode** curTreeLeaves = new CSGTreeNode*[pOctree->nMesh];
//        FacePair fPair;
//        ISectTriangle *curSurface;
//        TestTree testList;
//        std::vector<TMP_VInfo> points;
//
//        std::string randnumberout;
//        char str[64];
//
//        unsigned ibegin, iend;
//        switch (fDebugMode)
//        {
//        case -1:
//            ibegin = 0; iend = pOctree->nMesh;
//            break;
//        default:
//            ibegin = fDebugMode; iend = ibegin+1;
//            break;
//        }
//        for (unsigned i0 = ibegin; i0 < iend; i0++)
//        {
//            pMesh = pOctree->pMesh[i0];
//            int randid = rand()%pMesh->n_faces();
//            //randid = 31+offset++; // 11784, banana32
//            //randid = 392; //buuny
//            sprintf_s(str, "FLOOD:%u, rand:%d, faces:%d \n", i0, randid, pMesh->n_faces());
//            randnumberout += str;
//#ifdef DEBUG_FLOOD
//            randid = startId;
//#endif
//            curFace = pMesh->face_handle(randid);//46, sphere
//            
//            //AddTriangle(pMesh, curFace, GS::float4(1,0,0,1));
//            // 初始化第一个种子堆
//            seedQueueList.emplace();
//            SeedInfo &seedInfos = seedQueueList.back();
//            pMesh->property(pMesh->MarkPropHandle, curFace) = seeded;
//            seedInfos.relation = new Relation[pOctree->nMesh];
//            seedInfos.queue.emplace();
//            seedInfos.queue.back().seed = curFace;
//            seedInfos.queue.back().current = curFace;
//
//            auto surface = pMesh->property(pMesh->SurfacePropHandle, curFace);
//            memset(seedInfos.relation, 0, sizeof(Relation)*pOctree->nMesh);
//            seedInfos.relation[i0] = REL_SAME;
//            if (surface) MarkNARelation(surface, seedInfos.relation);
//
//            GetCorners(pMesh, curFace, v0, v1, v2);
//            Vec3d bc = (*v0+*v1+*v2)/3.0;
//            // 初始化独有的树
//            CSGTreeNode* curTree0 = copy2(pPosCSG->pRoot, curTreeLeaves);
//
//            for (unsigned i = 0; i < pOctree->nMesh; i++)
//            {
//                if (seedInfos.relation[i] == REL_UNKNOWN)
//                {
//                    bool IsInverse = pOctree->pMesh[i]->bInverse;
//                    if (!pOctree->pMesh[i]->BBox.Intersects(pMesh->BBox))
//                    {
//                        // must be out
//                        if (!IsInverse) curTreeLeaves[i]->relation = REL_OUTSIDE;
//                        else curTreeLeaves[i]->relation = REL_INSIDE;
//                        seedInfos.relation[i] = REL_NOT_AVAILABLE;
//                        continue;
//                    }
//
//                    seedInfos.relation[i] = PolyhedralInclusionTest(bc, pOctree, i, pOctree->pMesh[i]->bInverse);
//                    //sprintf_s(str, "PointInOut:%d, prim:%d, relation:%d \n", randid, i, seedInfos.relation[i]);
//                    //randnumberout += str;
//                }
//            }
//
//            auto rel0 = CompressCSGNodeIteration(curTree0);
//            if (rel0 != REL_NOT_AVAILABLE) continue;
//
//#ifdef DEBUG_FLOOD
//            int maxLimit = int(percent/100.0 * pMesh->n_faces());
//            int faceCount = 0;
//#endif
//            while (!seedQueueList.empty())
//            {
//                // 选取一个种子堆
//                auto &seeds = seedQueueList.front();
//
//                while (!seeds.queue.empty())
//                {
//                    while (!seeds.queue.empty())
//                    {
//                        if (pMesh->property(pMesh->MarkPropHandle, seeds.queue.front()[1]) != 2)
//                            break;
//                        seeds.queue.pop();
//                    }
//                    if (seeds.queue.empty()) break;
//
//                    // 初始化一个独立的种子
//                    curFace = seeds.queue.front()[1];
//                    relatedFace = seeds.queue.front()[0];
//                    faceQueue.push(curFace);
//                    seeds.queue.pop();
//                    GetRelationTable(pMesh, curFace, relatedFace, seeds.relation, pOctree->nMesh, pOctree, curRelationTable);
//
//                    GS::float4 color;
//                    color[3] = 1.0f;
//                    color[0] = rand() / (float)RAND_MAX;
//                    color[1] = rand() / (float)RAND_MAX;
//                    color[2] = 1.5f-color[0]-color[1];
//
//                    // 初始化该种子对应的种子堆
//                    seedQueueList.emplace();
//                    seedQueueList.back().relation = curRelationTable;
//
//                    // 生成关系树
//                    // TO-DO:通过检查可以减少这个树所需要的生成(ABC理论)
//                    memset(curTreeLeaves, 0, sizeof(CSGTreeNode*)*pOctree->nMesh);
//                    curTree = copy2(curTree0, curTreeLeaves);
//                    testList.clear();
//                    curRelation = ParsingCSGTree(pMesh, curRelationTable, pOctree->nMesh, curTree, curTreeLeaves, testList); // 未检查testList
//                    assert(!testList.size() || !testList.begin()->testTree->Parent);
//                    curSurface = pMesh->property(pMesh->SurfacePropHandle, curFace);
//#ifdef DEBUG_ISECT
//                    if (curSurface && (curSurface->segs.size() || curSurface->coplanarTris.size())) // 复合模式
//                    {
//                        color = GS::float4(0,0,0,1);
//                        for (auto &pair1:curSurface->segs)
//                        {
//                            switch (pair1.first)
//                            {
//                            case 0:
//                                color.x = 1;
//                                break;
//                            case 1:
//                                color.y = 1;
//                                break;
//                            case 2:
//                                color.z = 1;
//                                break;
//                            default:
//                                break;
//                            }
//                        }
//                        if (curSurface->segs.size() < 2)
//                        {
//                            for (auto& pair2: curSurface->coplanarTris)
//                            {
//                                switch (pair2.first)
//                                {
//                                case 0:
//                                    color.x = 1;
//                                    break;
//                                case 1:
//                                    color.y = 1;
//                                    break;
//                                case 2:
//                                    color.z = 1;
//                                    break;
//                                default:
//                                    break;
//                                }
//                            }
//                        }
//                    }
//                    color = GS::float4(1,0,0,1);
//#endif
//#ifdef FIXED_COLOR
//                    color = colorBuffer[icolor++];
//#endif
//                    if (curSurface && (curSurface->segs.size() || curSurface->coplanarTris.size())) // 复合模式
//                    {
//                        while (!faceQueue.empty())
//                        {                    
//                            while (!faceQueue.empty())
//                            {
//                                if (pMesh->property(pMesh->MarkPropHandle, faceQueue.front()) != 2)
//                                    break;
//                                faceQueue.pop();
//                            }
//                            if (faceQueue.empty()) break;
//                            curFace = faceQueue.front();
//#ifdef _DEBUG
//                            GetCorners(pMesh, curFace, v0, v1, v2);
//                            countd4 ++;
//#endif
//                            auto seedSurface = pMesh->property(pMesh->SurfacePropHandle, curFace);
//                            faceQueue.pop();
//
//                            // add neighbor
//                            bool IsSeed = false;
//                            ffItr = pMesh->ff_iter(curFace);
//                            int *markPtr;
//                            for (int i = 0; i < 3; i++, ffItr++)
//                            {
//                                if (!ffItr.is_valid()) break;
//                                markPtr = &(pMesh->property(pMesh->MarkPropHandle, *ffItr));
//                                if (*markPtr != 2)
//                                {
//                                    auto legSurface = pMesh->property(pMesh->SurfacePropHandle, *ffItr);
//                                    if (*markPtr != 1 && legSurface && CompareRelationSpace(seedSurface, legSurface))
//                                    {
//                                        *markPtr = 1; // queued
//                                        faceQueue.push(*ffItr);
//                                    }
//                                    else if (*markPtr != seeded)
//                                    {
//                                        FacePair fh;
//                                        fh[0] = curFace;
//                                        fh[1] = *ffItr;
//                                        seedQueueList.back().queue.push(fh);
//                                        *markPtr = seeded;
//                                        IsSeed = true;
//                                    }
//                                }
//                            }
//                            if (curRelation != REL_INSIDE || IsSeed)
//                            {
//                                points.clear();
//                                ParsingFace1(pMesh, curFace, pOctree->pMesh, points);
//                            }
//#ifdef DEBUG_ISECT
//                            AddTriangle(pMesh, curFace, color);
//#else
//#ifdef DEBUG_FLOOD
//                            faceCount ++;
//                            if (faceCount < maxLimit)
//                            {
//                                ParsingFace(pMesh, curFace, &testList, curRelation, pOctree->pMesh, points, pOctree, result, &color);
//                            }
//                            else
//                            {
//                                AddTriangle(pMesh, curFace, GS::float4(0.1,0.1,0.1,1));
//                            }
//#else
//                            ParsingFace(pMesh, curFace, &testList, curRelation, pOctree->pMesh, points, pOctree, result, &GS::float4(1.0f, 0.0f, 0.0f, 1.0f));
//#endif
//#endif
//                            pMesh->property(pMesh->MarkPropHandle, curFace) = 2; // processed
//                        }
//                    }
//                    else  // 简单模式
//                    {
//                        while (!faceQueue.empty())
//                        {            
//                            while (!faceQueue.empty())
//                            {
//                                if (pMesh->property(pMesh->MarkPropHandle, faceQueue.front()) != 2)
//                                    break;
//                                faceQueue.pop();
//                            }
//                            if (faceQueue.empty()) break;
//                            curFace = faceQueue.front();
//#ifdef _DEBUG
//                            GetCorners(pMesh, curFace, v0, v1, v2);
//                            countd3 ++;
//                            assert(curRelation == REL_INSIDE || curRelation == REL_SAME);
//#endif
//
//                            faceQueue.pop();
//#ifdef DEBUG_ISECT
//                            AddTriangle(pMesh, curFace);
//#else
//#ifdef DEBUG_FLOOD
//                            faceCount ++;
//                            if (faceCount < maxLimit)
//                            {
//                                if (curRelation == REL_SAME) AddTriangle(pMesh, curFace, color);
//                            }
//                            else
//                            {
//                                AddTriangle(pMesh, curFace, GS::float4(0.1,0.1,0.1,1));
//                            }
//#else
//                            if (curRelation == REL_SAME) AddTriangle(pMesh, curFace);//AddTriangle(pMesh, curFace, color);
//#endif
//#endif
//                            pMesh->property(pMesh->MarkPropHandle, curFace) = 2; // processed
//
//                            // add neighbor
//                            ffItr = pMesh->ff_iter(curFace);
//                            int *markPtr;
//                            for (int i = 0; i < 3; i++, ffItr++)
//                            {
//                                if (!ffItr.is_valid()) break;
//                                markPtr = &(pMesh->property(pMesh->MarkPropHandle, *ffItr));
//                                if (*markPtr != 2)
//                                {
//                                    auto legSurface = pMesh->property(pMesh->SurfacePropHandle, *ffItr);
//                                    if (*markPtr != 1 && (!legSurface || (!legSurface->segs.size() && !legSurface->coplanarTris.size())))
//                                    {
//                                        *markPtr = 1; // queued
//                                        faceQueue.push(*ffItr);
//                                    }
//                                    else if (*markPtr != seeded)
//                                    {
//                                        FacePair fh;
//                                        fh[0] = curFace;
//                                        fh[1] = *ffItr;
//                                        seedQueueList.back().queue.push(fh);
//                                        *markPtr = seeded;
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//                seedQueueList.pop();
//            }
//            delete curTree0;
//        }
//        //StdOutput(randnumberout.c_str());
//        delete [] curTreeLeaves;
//    }
}

