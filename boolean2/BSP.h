#pragma once
#include <vector>
#include "preps.h"
#include "xmemory.h"

namespace Boolean
{
    class BSPTree
    {
    public:
        struct Node
        {
            XPlane sp;
            Oriented_side relation;
            Node *posChind = nullptr, *negChild = nullptr;

            //std::vector<IPolygon*> coins;

            bool isLeaf() const { return posChind == nullptr; }
        };

    public:
        BSPTree() {}
        ~BSPTree() { reset(); }

        void buildNoCross(std::vector<Triangle*>& faces);
        //void build(std::vector<IPolygon*>& faces);
        void reset();
        Oriented_side classify(const MyVertex& v, XPlane* bspPlane = 0) const;

    private:
        Node* buildRecursion(std::vector<Triangle*>& faces);

        Node *mp_root = nullptr;

    };

}
