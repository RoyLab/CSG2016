#pragma once
#include <vector>
#include "preps.h"
#include "RegularMesh.h"

namespace Boolean
{
	class MemoryManager
	{
		typedef cyPointT VPoint;
	public:
		std::vector<XPlaneBase> planes;
		std::vector<MyEdge> edges;
		std::vector<MyVertex> vertices;

		std::vector<VPoint>	points;
		std::vector<XPoint>	ppoints;

	public:
		~MemoryManager() {}
		static MemoryManager* getInstance();

		int insertVertices(VPoint* begin, VPoint* end);
		uint32_t getEdgeId(uint32_t a, uint32_t b, IPolygon* facePtr);

		//// access
		//XPlane* getPlaneBuffer() { return m_pplanes; }
		//const XPlane* getPlaneBuffer() const { return m_pplanes; }

	private:
		MemoryManager() {}
	};

	inline const MyEdge& xcedge(uint32_t id) { return MemoryManager::getInstance()->edges[id]; }
	inline MyEdge& xedge(uint32_t id) { return MemoryManager::getInstance()->edges[id]; }
	inline const std::vector<MyEdge>& xedges() { return MemoryManager::getInstance()->edges; }

	inline const MyVertex& xcvertex(uint32_t id) { return MemoryManager::getInstance()->vertices[id]; }
	inline MyVertex& xvertex(uint32_t id) { return MemoryManager::getInstance()->vertices[id]; }
	inline const std::vector<MyVertex>& xvertices() { return MemoryManager::getInstance()->vertices; }

	inline const XPlaneBase& xplane(uint32_t id) { return MemoryManager::getInstance()->planes[id]; }
}
