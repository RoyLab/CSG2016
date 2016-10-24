#include "precompile.h"
#include "preps.h"
#include "xmemory.h"

namespace Boolean
{
	const XPlaneBase & XPlane::base() const
	{
		return MemoryManager::getInstance()->planes.at(std::abs(id) - 1);
	}
}