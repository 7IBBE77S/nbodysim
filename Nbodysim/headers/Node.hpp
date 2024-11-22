#pragma once
#include "Vec2.hpp"
#include "Quad.hpp"

struct Node
{

    // struct Range
    // {
    //     size_t start;
    //     size_t end;
    //     size_t size() const { return end - start; }
    // };
    size_t children;
    size_t next;
    Vec2 pos;
    float mass;
    Quad quad;
    size_t bodies_start;
    size_t bodies_end;
    size_t body_count;

    Node(size_t next, Quad quad)
        : children(0), next(next), pos(Vec2::zero()), mass(0.0f), quad(quad), bodies_start(0), bodies_end(0), body_count(0) {} // Initialize body_count

    bool is_leaf() const { return children == 0; }
    bool is_branch() const { return children != 0; }
    bool is_empty() const { return mass == 0.0f; }
};