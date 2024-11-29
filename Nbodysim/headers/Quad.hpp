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
    Quad(const Quad& other) : center(other.center), size(other.size) {}
    
    Quad(Quad&& other) noexcept : center(std::move(other.center)), size(other.size) {}

       Quad& operator=(const Quad& other) {
        center = other.center;
        size = other.size;
        return *this;
    }

    Quad& operator=(Quad&& other) noexcept {
        center = std::move(other.center);
        size = other.size;
        return *this;
    }


    static Quad new_containing(const std::vector<Body>& bodies) {
        Vec2 min_bounds = Vec2::broadcast(std::numeric_limits<float>::max());
        Vec2 max_bounds = Vec2::broadcast(std::numeric_limits<float>::lowest());

        for (const auto& body : bodies) {
            min_bounds = min_bounds.min_by_component(body.pos);
            max_bounds = max_bounds.max_by_component(body.pos);
        }

        Vec2 center = (min_bounds + max_bounds) * 0.5f;
        Vec2 size_vec = max_bounds - min_bounds;
        float size = size_vec.component_max();

        return Quad(center, size);
    }

    size_t find_quadrant(Vec2 pos) const {
        return ((pos.y > center.y) << 1) | (pos.x > center.x);
    }

    Quad into_quadrant(size_t quadrant) const {
        float new_size = size * 0.5f;
        Vec2 offset = Vec2::unit_x() * ((quadrant & 1) - 0.5f) + 
                     Vec2::unit_y() * ((quadrant >> 1) - 0.5f);
        Vec2 new_center = center + offset * new_size;
        return Quad(new_center, new_size);
    }

    std::array<Quad, 4> subdivide() const {
        return {into_quadrant(0), into_quadrant(1), into_quadrant(2), into_quadrant(3)};
    }
};