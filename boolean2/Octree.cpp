#include "Octree.h"
#include "CGALext.h"
#include "RegularMesh.h"


namespace Boolean
{
    int Octree::MAX_TRIANGLE_COUNT = 30;
    int Octree::MAX_LEVEL = 9;

    void Octree::release()
    {
        SAFE_DELETE(mp_root);
        mp_meshes = nullptr;
        m_nMesh = 0;
    }

    Octree::Node* Octree::createRootNode()
    {
        Node* root = new Node;
        CGAL::Bbox_3 tmpBox = mp_meshes[0]->get_bbox();

        for (size_t i = 0; i < m_nMesh; i++)
        {
            auto pcMesh = mp_meshes[i];
            tmpBox += pcMesh->get_bbox();

            for (auto fItr = pcMesh->facets_begin(); fItr != pcMesh->facets_end(); fItr++)
            {
                root->triCount++;
                auto &pTable = root->triTable[i];
                if (!pTable)
                    pTable = new TriList;

                pTable->push_back(fItr);
            }
        }

        root->bbox = enlarge(tmpBox, 1e-5);
        return root;
    }

    void Octree::build(const std::vector<RegularMesh*>& meshList, std::vector<Node*>* isectNodes)
    {
        if (!meshList.size()) throw std::exception("mesh list size 0.");

        mp_meshes = meshList.data();
        m_nMesh = meshList.size();

        mp_root = createRootNode();
        build(mp_root, 0, isectNodes);
    }

    void Octree::build(Node* root, size_t level, std::vector<Node*>* isectNodes)
    {
        assert(root);

        if (root->triTable.size() <= 1)
        {
            root->type = NODE_SIMPLE;
        }
        else if (root->triCount <= MAX_TRIANGLE_COUNT || level > MAX_LEVEL)
        {
            root->type = NODE_COMPOUND;
            if (isectNodes) isectNodes->push_back(root);
        }
        else
        {
            root->type = NODE_MIDSIDE;
            root->pChildren = new Node[8];

            Vector minOffset, maxOffset;
            Bbox bbox = root->bbox;
            Vector step = (bbox.max() - bbox.min()) * 0.5;
    
            for (unsigned i = 0; i < 8 ; i++)
            {
                auto pChild = &root->pChildren[i];

                maxOffset = Vector(i & 4 ? 0 : -step[0],
                    i & 2 ? 0 : -step[1],
                    i & 1 ? 0 : -step[2]);

                minOffset = Vector(i & 4 ? step[0] : 0,
                    i & 2 ? step[1] : 0,
                    i & 1 ? step[2] : 0);

                pChild->bbox = Bbox(bbox.min() + minOffset, 
                    bbox.max() + maxOffset);

                pChild->pParent = root;
            }

            for (auto &triTab: root->triTable)
            {
                size_t meshId = triTab.first;
                RegularMesh* pMesh = mp_meshes[meshId];
                TriList &parentMeshes = *triTab.second;

                const size_t tn = parentMeshes.size();
                for (size_t i = 0; i < tn; i++)
                {
                    auto fh = parentMeshes[i];
                    for (unsigned j = 0; j < 8; j++)
                    {
                        Node &child = root->pChildren[j];
                        TriList* triList = nullptr;
                        int isInbox = -1; // -1 no intersection, 1 bbox in box, 0 bbox not in box

                        if (is_inside_box(child.bbox, fh->data->bbox))
                            isInbox = 1;
                        else if (CGAL::do_intersect(child.bbox.bbox(), fh->data->triangle))
                            isInbox = 0;

                        if (isInbox >= 0)
                        {
                            if (!triList)
                                triList = child.triTable.emplace(meshId, new TriList()).first->second;

                            triList->push_back(fh);
                            child.triCount++;

                            if (isInbox == 1)
                                break;
                        }
                    }
                }
            }

            root->triTable.clear();
            for (size_t i = 0; i < 8; i++)
                build(root->pChildren + i, level + 1, isectNodes);

        }
    }


} // namespace CSG

