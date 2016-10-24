#pragma once
#include <vector>
#include "preps.h"
#include "RegularMesh.h"

namespace Boolean
{
	class MemoryManager
	{
		typedef cyPointT Point;
	public:
		std::vector<XPlaneBase> planes;
		std::vector<MyEdge> edges;
		std::vector<MyVertex> vertices;
		std::vector<Point>	points;

	public:
		~MemoryManager();
		static MemoryManager* getInstance();

		int insertVertices(Point* begin, Point* end);
		int getEdgeId(int a, int b, IPolygon* facePtr);

		//// access
		//XPlane* getPlaneBuffer() { return m_pplanes; }
		//const XPlane* getPlaneBuffer() const { return m_pplanes; }

	private:
		MemoryManager() {}
	};
}
