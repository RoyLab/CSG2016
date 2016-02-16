#pragma once
#include "csg.h"
#include "csgdefs.h"
#include "MPMesh.h"
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

    struct DiffMeshInfo
    {
        unsigned ID;
        Relation Rela;

        DiffMeshInfo(
            unsigned i, 
            Relation rel = REL_UNKNOWN
            ): ID(i), Rela(rel)
        {}
    };

    typedef std::map<size_t, std::vector<MPMesh::FaceHandle>> TriTableT;

    class Octree
    {
    public:
        struct Node
        {
            Bbox_3 bbox;
            NodeType type;

            Node *pChildren, *pParent;

            //std::vector<DiffMeshInfo> DiffMeshIndex;
            TriTableT   triTable;
            size_t      triCount;

            //void *pRelationData;

            Node();
            ~Node();
        };

    public:
        Octree(){}
        ~Octree(){ release(); }

        void build(const std::vector<MPMesh*>& meshList, std::vector<Node*>* leaves = nullptr);
        void release();

    private:
        Node*               mp_root;
        MPMesh* const*      mp_meshes;
        unsigned            m_nMesh;
    };


    //typedef Octree<MPMesh> MyOctree;
    //Octree<>* BuildOctree(MPMesh** meshList, unsigned nMesh);
    //Relation PolyhedralInclusionTest(Vec3d& point, Octree<>* pOctree, unsigned meshId, bool = false);

    //inline bool IsLeaf(OctreeNode* node) {return !node->Child;}

}// namespace CSG

