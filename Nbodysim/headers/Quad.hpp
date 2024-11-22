#pragma once
#include "Vec2.hpp"
#include "Body.hpp"
#include <vector>
#include <algorithm>
#include <array>

struct Quad {
    Vec2 center;
    float size;

    Quad() = default;
    Quad(Vec2 center, float size) : center(center), size(size) {}

    static Quad new_containing(const std::vector<Body>& bodies) {
        float min_x = std::numeric_limits<float>::max();
        float min_y = std::numeric_limits<float>::max();
        float max_x = std::numeric_limits<float>::lowest();
        float max_y = std::numeric_limits<float>::lowest();

        for (const auto& body : bodies) {
            min_x = std::min(min_x, body.pos.x);
            min_y = std::min(min_y, body.pos.y);
            max_x = std::max(max_x, body.pos.x);
            max_y = std::max(max_y, body.pos.y);
        }

        Vec2 center((min_x + max_x) * 0.5f, (min_y + max_y) * 0.5f);
        float size = std::max(max_x - min_x, max_y - min_y);

        return Quad(center, size);
    }

    size_t find_quadrant(Vec2 pos) const {
        return ((pos.y > center.y) << 1) | (pos.x > center.x);
    }

    Quad into_quadrant(size_t quadrant) const {
        float new_size = size * 0.5f;
        Vec2 new_center = center;
        new_center.x += ((quadrant & 1) - 0.5f) * new_size;
        new_center.y += ((quadrant >> 1) - 0.5f) * new_size;
        return Quad(new_center, new_size);
    }

    std::array<Quad, 4> subdivide() const {
        return {into_quadrant(0), into_quadrant(1), into_quadrant(2), into_quadrant(3)};
    }
};