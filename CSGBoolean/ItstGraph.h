#pragma once
#include <map>
#include <deque>

#include "MyMesh.h"
#include "macroutil.h"
#include "AlgUserData.h"
#include "ItstAlg.h"

namespace CSG
{
    /* 最小可运行版本，只考虑单折痕，不考虑孤立顶点 */
    class ItstGraph
    {
        DEFINE_HANDLES;
        enum {SEQ = 0x01, INV = 0x02};

        struct Edge;
        struct Adj
        {
            int vId = -1;
            int eId = -1;
            bool visited = false;
        };

        struct Node
        {
            VProxyItr vproxy;
            std::deque<Adj> edges;
        };

        struct Edge
        {
            int vId[2];
            int direction = 0;
        };

        typedef std::map<int, int> NodeMap;

        COMMON_PROPERTY(bool, bValid);
        COMMON_PROPERTY(std::deque<Node>, nodes);
        COMMON_PROPERTY(std::deque<Edge>, edges);
        COMMON_PROPERTY(NodeMap, maps);
        COMMON_PROPERTY_POINTER(ItstAlg, alg);
    public:

        ItstGraph(FH fh, ItstAlg* data, int meshId);
        ~ItstGraph();

        void addEdge(Edge& e);
        void addNode(Node& node);

    private:
        void sortEdgeDirection();
    };


}