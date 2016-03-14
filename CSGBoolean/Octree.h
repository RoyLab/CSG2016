#pragma once

#include <vector>
#include <map>

#include "csg.h"
#include "csgdefs.h"
#include "macroutil.h"
#include "MyMesh.h"



namespace CSG
{
    enum NodeType
    {
        NODE_UNKNOWN = 0,
        NODE_MIDSIDE,
        NODE_SIMPLE,
        NODE_COMPOUND
    };

    class Octree
    {
        typedef CGAL::Point_3<K> Point;
        typedef CGAL::Vector_3<K> Vector;
        typedef Cube_3 Bbox;
        typedef std::vector<MyMesh::Face_handle> TriList;
        typedef std::map<uint32_t, TriList*> TriTableT;

    public:
        struct Node
        {
            Bbox        bbox;
            NodeType    type;

            Node        *pChildren = nullptr,
                        *pParent = nullptr;

            TriTableT   triTable;
            size_t      triCount;

            virtual ~Node() { 
                SAFE_DELETE_ARRAY(pChildren); 
                for (auto pair : triTable)
                    SAFE_DELETE(pair.second);
            }
        };

    public:
        Octree(){}
        ~Octree(){ release(); }

        void build(const std::vector<MyMesh*>& meshList, std::vector<Node*>* isectNodes = nullptr);
        void release(); // not release meshlist
        Node* getRoot() const { return mp_root; }

    private:
        Node* createRootNode();
        void build(Node* root, size_t level, std::vector<Node*>* isectNodes = nullptr);


    private:
        Node*               mp_root;
        MyMesh* const*      mp_meshes;
        size_t              m_nMesh;

        STATIC_PROPERTY(int, MAX_TRIANGLE_COUNT);
        STATIC_PROPERTY(int, MAX_LEVEL);
    };

}// namespace CSG

