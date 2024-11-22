#pragma once
#include "Vec2.hpp"

struct Body {
    Vec2 pos;
    Vec2 vel;
    Vec2 acc;
    float mass;
    float radius;

    Body() = default;
    Body(Vec2 pos, Vec2 vel, float mass, float radius)
        : pos(pos), vel(vel), acc(Vec2::zero()), mass(mass), radius(radius) {}

    void update(float dt) {
        vel += acc * dt;
        pos += vel * dt;
    }
};