#pragma once
#include "csg.h"
#include "csgdefs.h"
#include "macro_util.h"
#include "MyMesh.h"

#include <vector>
#include <map>


namespace CSG
{
    enum NodeType
    {
        NODE_UNKNOWN = 0,
        NODE_MIDSIDE,
        NODE_SIMPLE,
        NODE_COMPOUND
    };

    //struct DiffMeshInfo
    //{
    //    unsigned ID;
    //    Relation Rela;

    //    DiffMeshInfo(
    //        unsigned i, 
    //        Relation rel = REL_UNKNOWN
    //        ): ID(i), Rela(rel)
    //    {}
    //};


    class Octree
    {
        typedef CGAL::Point_3<K> Point;
        typedef CGAL::Vector_3<K> Vector;
        typedef Cube_3 Bbox;
    public:
        struct Node
        {
            typedef std::map<size_t, std::vector<MyMesh::Face_handle>> TriTableT;

            Bbox        bbox;
            NodeType    type;

            Node        *pChildren = nullptr,
                        *pParent = nullptr;

            TriTableT   triTable;
            size_t      triCount;

            virtual ~Node() { SAFE_DELETE_ARRAY(pChildren); }

            //std::vector<DiffMeshInfo> DiffMeshIndex;
            //void *pRelationData;
        };

    public:
        Octree(){}
        ~Octree(){ release(); }

        void build(const std::vector<MyMesh*>& meshList, std::vector<Node*>* leaves = nullptr);
        void release();

    private:
        Node* createRootNode();
        void build(Node* root, size_t level);

    private:
        Node*               mp_root;
        MyMesh* const*      mp_meshes;
        unsigned            m_nMesh;

        static int MAX_TRIANGLE_COUNT;
        static int MAX_LEVEL;
    };


    //typedef Octree<MyMesh> MyOctree;
    //Octree<>* BuildOctree(MyMesh** meshList, unsigned nMesh);
    //Relation PolyhedralInclusionTest(Vec3d& point, Octree<>* pOctree, unsigned meshId, bool = false);

    //inline bool IsLeaf(OctreeNode* node) {return !node->Child;}

}// namespace CSG

