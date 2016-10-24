#pragma once

#include <vector>
#include <map>

#include <macros.h>
#include "global.h"
#include "csg.h"
#include "csgdefs.h"



namespace Boolean
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
        typedef CGAL::Point_3<Depick> Point;
        typedef CGAL::Vector_3<Depick> Vector;
        typedef CGAL::Iso_cuboid_3<Depick> NodeShape;
        typedef std::vector<size_t> TriList;
        typedef std::map<uint32_t, TriList*> TriTableT;

    public:
        struct Node
        {
            NodeShape   bbox;
            NodeType    type = NODE_UNKNOWN;

            Node        *pChildren = nullptr,
                        *pParent = nullptr;

            TriTableT   triTable;
            size_t      triCount = 0;

            virtual ~Node() { 
                SAFE_DELETE_ARRAY(pChildren); 
                for (auto pair : triTable)
                    SAFE_DELETE(pair.second);
            }
        };

    public:
        Octree(){}
        ~Octree(){ release(); }

        void build(const std::vector<RegularMesh*>& meshList, std::vector<Node*>* isectNodes = nullptr);
        void release(); // not release meshlist
        Node* getRoot() const { return mp_root; }

    private:
        Node* createRootNode();
        void build(Node* root, size_t level, std::vector<Node*>* isectNodes = nullptr);


    private:
        Node*					mp_root = nullptr;
		RegularMesh* const*     mp_meshes = nullptr;
        size_t					m_nMesh = 0;

        STATIC_PROPERTY(int, MAX_TRIANGLE_COUNT);
        STATIC_PROPERTY(int, MAX_LEVEL);
    };

}// namespace CSG

