#pragma once
#include <map>
#include <deque>

#include "MyMesh.h"
#include "macroutil.h"
#include "AlgUserData.h"
#include "ItstAlg.h"
#include "csgdefs.h"
#include "plane_reps.h"

namespace CSG
{
    /* 最小可运行版本，只考虑单折痕，不考虑孤立顶点 */
    class ItstGraph
    {
        DEFINE_HANDLES;
        enum {SEQ = 0x01, INV = 0x02};

    public:
        struct Edge;
        struct Adj
        {
            int vId = -1; // local id
            int eId = -1; // local id
            bool visited = false;

            K::Vector_3 calcDirection(ItstGraph* graph, const PBPoint<K>& orgin);
        };

        struct Node
        {
            VProxyItr vproxy;
            int id = -1; // local id
            std::deque<Adj> edges;
            IIndicatorVector* indicator = nullptr;

            int mark = UNVISITED;

            ~Node() { SAFE_DELETE(indicator); }
        };

        struct Edge
        {
            int vId[2]; // global id
            int direction = 0;
        };

        typedef std::map<int, int> NodeMap;
        typedef std::vector<int> IdContainer;

        COMMON_PROPERTY(bool, bValid);
        COMMON_PROPERTY(std::deque<Node>, nodes);
        COMMON_PROPERTY(std::deque<Edge>, edges);
        COMMON_PROPERTY(NodeMap, maps);
        COMMON_PROPERTY(IdContainer, ids);
        COMMON_PROPERTY_POINTER(ItstAlg, alg);
    public:

        ItstGraph(FH fh, ItstAlg* data, int meshId);
        ~ItstGraph();

        template <class Container>
        void floodFilling(VH vh, SampleIndicatorVector& sample, Container& ids)
        {
            for (int &id : ids)
                m_ids.push_back(id);

            auto result = m_maps.find(vh->data->proxy->pointer()->idx);
            assert(result != m_maps.end());

            int startId = result->second;
            m_nodes[startId].indicator = new SampleIndicatorVector(sample);

            floodFilling(startId);
        }

        void floodFilling(int startId);

        void addEdge(Edge& e);
        void addNode(Node& node);
        bool isNodeVisited(Adj& adj) const { m_nodes[adj.vId].mark == VISITED; }

    private:
        void sortEdgeDirection();
    };


}