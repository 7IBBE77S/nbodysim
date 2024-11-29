#pragma once
#include "Vec2.hpp"
#include <algorithm>


struct alignas(16) AABB
{
    Vec2 min;
    Vec2 max;


    AABB combine(const AABB &other) const
    {
        return AABB{
            Vec2(std::min(min.x, other.min.x), std::min(min.y, other.min.y)),
            Vec2(std::max(max.x, other.max.x), std::max(max.y, other.max.y))};
    }
    bool intersects(const AABB &other) const
    {
        return (this->min.x <= other.max.x && this->max.x >= other.min.x) &&
               (this->min.y <= other.max.y && this->max.y >= other.min.y);
    }
    static AABB from_body(const Body &body)
    {
        Vec2 r(body.radius, body.radius);
        return {body.pos - r, body.pos + r};
    }

  
};