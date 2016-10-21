#pragma once
#include "global.h"

namespace Boolean
{
	class MemoryManager
	{
	public:
		~MemoryManager();
		static MemoryManager* getInstance();

		// access
		XPlane* getPlaneBuffer() { return m_planes; }
		const XPlane* getPlaneBuffer() const { return m_planes; }

	private:
		MemoryManager() {}

		XPlane *m_planes;
		XPoint *m_points;
	};
}
