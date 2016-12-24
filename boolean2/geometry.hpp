#pragma once

#include "csgdefs.h"

namespace Boolean
{
    inline cyPointT get_projection_vector_unit(const cyPointT& p)
    {
        int max_component = 0;
        Real max_value = std::abs(p.x);
        
        Real tmp = std::abs(p.y);
        if (tmp > max_value)
        {
            max_value = tmp;
            max_component = 1;
        }

        tmp = std::abs(p.z);
        if (tmp > max_value)
        {
            max_value = tmp;
            max_component = 2;
        }

        cyPointT result(0, 0, 0);
        result[max_component] = std::copysign(1, p[max_component]);
        return result;
    }
}

