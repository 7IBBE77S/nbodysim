#pragma once
#include "Vec2.hpp"
#include <algorithm>

struct AABB {
    Vec2 min;
    Vec2 max;
    

    bool overlaps(const AABB& other) const {
        return (max.x > other.min.x && min.x < other.max.x) &&
               (max.y > other.min.y && min.y < other.max.y);
    }

    static AABB combine(const AABB& a, const AABB& b) {
        return {
            {std::min(a.min.x, b.min.x), std::min(a.min.y, b.min.y)},
            {std::max(a.max.x, b.max.x), std::max(a.max.y, b.max.y)}
        };
    }
};