#pragma once
#include "Vec2.hpp"
#include "Quad.hpp"
#include <vector>
#include <algorithm>
#include <array>
#include <iostream>


struct Range
{
    size_t start;
    size_t end;

    Range(size_t start = 0, size_t end = 0) : start(start), end(end) {}
    Range(const Range &other) : start(other.start), end(other.end) {}

    Range &operator=(const Range &other)
    {
        if (this != &other)
        {
            start = other.start;
            end = other.end;
        }
        return *this;
    }

    size_t size() const { return end - start; }
};

struct alignas(32) Node
{
    // put frequently accessed data grouped together for cache efficiency
    // Hot data together
    struct alignas(16)
    {
        Vec2 pos;   // 8 bytes
        float mass; // 4 bytes
        Quad quad;  // 12 bytes
    } data;          // total 24 byte aligned
                     // Cold data
    size_t children; // 8 bytes - Index to first child
    size_t next;     // 8 bytes - Index to next node
    Range bodies;    // 16 bytes - Body range using semantic struct
    size_t depth;

    Node(size_t next = 0, Quad quad = Quad(), size_t depth = 0)
        : data{Vec2::zero(), 0.0f, quad}, children(0), next(next), bodies{0, 0}, depth(depth) {}

    bool is_leaf() const { return children == 0; }
    bool is_branch() const { return children != 0; }
    bool is_empty() const { return data.mass == 0.0f; }
};

// class NodeSoA
// {
// private:
//     // Hot data cluster
//     struct HotData
//     {
//         std::vector<Vec2> positions;
//         std::vector<float> masses;
//         std::vector<Quad> quads;
//     };

//     // Cold data cluster
//     struct TreeData
//     {
//         std::vector<size_t> children;
//         std::vector<size_t> next_indices;
//         std::vector<Range> body_ranges;
//     };

//     HotData hot;
//     TreeData tree;

//     size_t size_;

// public:
//     NodeSoA(size_t capacity = 1024) : size_(0)
//     {
//         reserve(capacity);
//     }
//     void set_bodies(size_t idx, const Range &range)
//     {
//         tree.body_ranges[idx] = range;
//     }

//     void reserve(size_t capacity)
//     {
//         hot.positions.reserve(capacity);
//         hot.masses.reserve(capacity);
//         hot.quads.reserve(capacity);
//         tree.children.reserve(capacity);
//         tree.next_indices.reserve(capacity);
//         tree.body_ranges.reserve(capacity);
//     }
//     void clear()
//     {
//         hot.positions.clear();
//         hot.masses.clear();
//         hot.quads.clear();
//         tree.children.clear();
//         tree.next_indices.clear();
//         tree.body_ranges.clear();
//         size_ = 0;
//     }

//     size_t add_node(size_t next = 0, Quad quad = Quad())
//     {
//         if (size_ >= hot.positions.capacity())
//         {
//             size_t new_capacity = (hot.positions.capacity() == 0) ? 1 : hot.positions.capacity() * 2;
//             reserve(new_capacity);
//         }

//         size_t idx = size_++;
//         hot.positions.push_back(Vec2::zero());
//         hot.masses.push_back(0.0f);
//         hot.quads.push_back(quad);
//         tree.children.push_back(0);
//         tree.next_indices.push_back(next);
//         tree.body_ranges.push_back(Range{0, 0});
//         return idx;
//     }

//     bool is_leaf(size_t idx) const { return tree.children[idx] == 0; }
//     bool is_branch(size_t idx) const { return tree.children[idx] != 0; }
//     bool is_empty(size_t idx) const { return hot.masses[idx] == 0.0f; }
//     bool is_empty() const { return size_ == 0; }

//     const Vec2 &position(size_t idx) const { return hot.positions[idx]; }
//     const float &mass(size_t idx) const { return hot.masses[idx]; }
//     const Quad &quad(size_t idx) const { return hot.quads[idx]; }
//     const size_t &child(size_t idx) const { return tree.children[idx]; }
//     const size_t &next(size_t idx) const { return tree.next_indices[idx]; }
//     const Range &bodies(size_t idx) const { return tree.body_ranges[idx]; }
//     const size_t &depth(size_t idx) const { return tree.children[idx]; }

//     Vec2 &position(size_t idx) { return hot.positions[idx]; }
//     float &mass(size_t idx) { return hot.masses[idx]; }
//     Quad &quad(size_t idx) { return hot.quads[idx]; }
//     size_t &child(size_t idx) { return tree.children[idx]; }
//     size_t &next(size_t idx) { return tree.next_indices[idx]; }
//     Range &bodies(size_t idx) { return tree.body_ranges[idx]; }
//     size_t &depth(size_t idx) { return tree.children[idx]; }

//     size_t size() const { return size_; }
// };