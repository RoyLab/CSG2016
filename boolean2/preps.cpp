#include "precompile.h"
#include "preps.h"
#include "xmemory.h"
#include "adaptive.h"

namespace Boolean
{
	const XPlaneBase & XPlane::base() const
	{
		return xplane(std::abs(id) - 1);
	}

	void XLine::makePositive(XPlane & input)
	{
		const Real* mat[3] = { m_planes[0].data(), m_planes[1].data(), input.data() };
		if (adaptiveDet3x3Sign(mat) < 0.0)
			input.inverse();
	}
}