#pragma once
#include "Body.hpp"
#include "Quadtree.hpp"
#include "AABB.hpp"
#include <vector>
#include <random>
#include <numeric>
#include <cmath>
#include <numbers> 

struct BVHNode
{
    AABB bounds;
    BVHNode *left = nullptr;
    BVHNode *right = nullptr;
    size_t bodyIndex = SIZE_MAX;

    bool isLeaf() const { return left == nullptr && right == nullptr; }
    ~BVHNode()
    {
        delete left;
        delete right;
    }
};

// struct alignas(16) Body {
//     Vec2 pos;     // 8 bytes
//     Vec2 vel;     // 8 bytes
//     Vec2 acc;     // 8 bytes
//     float mass;   // 4 bytes
//     float radius; // 4 bytes
// };

// // Add before Simulation class
// struct BodySystem {
//     static constexpr size_t SIMD_WIDTH = 4;
//     std::vector<Vec2> positions;    // SOA layout
//     std::vector<Vec2> velocities;
//     std::vector<Vec2> accelerations;
//     std::vector<float> masses;
//     std::vector<float> radii;

//     void resize(size_t n) {
//         positions.resize(n);
//         velocities.resize(n);
//         accelerations.resize(n);
//         masses.resize(n);
//         radii.resize(n);
//     }

//     // Convert from AOS to SOA
//     void fromBodies(const std::vector<Body>& bodies) {
//         resize(bodies.size());
//         for(size_t i = 0; i < bodies.size(); i++) {
//             positions[i] = bodies[i].pos;
//             velocities[i] = bodies[i].vel;
//             accelerations[i] = bodies[i].acc;
//             masses[i] = bodies[i].mass;
//             radii[i] = bodies[i].radius;
//         }
//     }
// };

// inline uint32_t calculateMortonCode(Vec2 pos) {
//     // Normalize positions to [0,1] range
//     // Assuming reasonable bounds for your simulation
//     const float scale = 1.0f / 1024.0f;  // Adjust based on your world size
//     uint32_t x = static_cast<uint32_t>(pos.x * scale) & 0x0000FFFF;
//     uint32_t y = static_cast<uint32_t>(pos.y * scale) & 0x0000FFFF;

//     // Interleave bits of x and y to create Morton code
//     x = (x | (x << 8)) & 0x00FF00FF;
//     x = (x | (x << 4)) & 0x0F0F0F0F;
//     x = (x | (x << 2)) & 0x33333333;
//     x = (x | (x << 1)) & 0x55555555;

//     y = (y | (y << 8)) & 0x00FF00FF;
//     y = (y | (y << 4)) & 0x0F0F0F0F;
//     y = (y | (y << 2)) & 0x33333333;
//     y = (y | (y << 1)) & 0x55555555;

//     return x | (y << 1);
// }
// // Helper function to calculate surface area of an AABB
// inline float surfaceArea(const AABB& bounds) {
//     Vec2 diff = bounds.max - bounds.min;
//     return diff.x * diff.y;  // For 2D, this is actually the area
// }

class Simulation
{
public:
    float dt;
    size_t frame;
    std::vector<Body> bodies;
    Quadtree quadtree;
    std::vector<size_t> bodyIndices;

    Simulation()
        : dt(0.2f), frame(0), quadtree(1.0f, 1.0f, 16)
    {
        size_t n = 8000;
        bodies.reserve(n);
        bodyIndices.reserve(n);
        bodies = uniform_disc(n);
    }

    void step()
    {
        iterate();
        collide();
        attract();
        ++frame;
    }

private:
    void iterate()
    {
        for (auto &body : bodies)
        {
            body.update(dt);
        }
    }

    void attract()
    {
        quadtree.build(bodies);
        for (auto &body : bodies)
        {
            body.acc = quadtree.acc(body.pos);
        }
    }

    void collide()
    {
        if (bodies.empty())
            return;

        bodyIndices.resize(bodies.size());
        std::iota(bodyIndices.begin(), bodyIndices.end(), 0);

        BVHNode *root = buildBVH(bodyIndices, 0, bodies.size());
        collideBVH(root, root);
        delete root;
    }

    AABB computeBodyAABB(const Body &body)
    {
        Vec2 r(body.radius, body.radius);
        return {body.pos - r, body.pos + r};
    }

    BVHNode *buildBVH(std::vector<size_t> &indices, size_t start, size_t end)
    {
        if (start >= end)
            return nullptr;

        BVHNode *node = new BVHNode();

        if (end - start == 1)
        {
            node->bodyIndex = indices[start];
            node->bounds = computeBodyAABB(bodies[node->bodyIndex]);
            return node;
        }

        AABB bounds = computeBodyAABB(bodies[indices[start]]);
        for (size_t i = start + 1; i < end; ++i)
        {
            bounds = AABB::combine(bounds, computeBodyAABB(bodies[indices[i]]));
        }
        node->bounds = bounds;

        size_t mid = start + (end - start) / 2;
        bool splitX = (bounds.max.x - bounds.min.x) > (bounds.max.y - bounds.min.y);

        std::sort(indices.begin() + start, indices.begin() + end,
                  [&](size_t a, size_t b)
                  {
                      return splitX ? bodies[a].pos.x < bodies[b].pos.x : bodies[a].pos.y < bodies[b].pos.y;
                  });

        node->left = buildBVH(indices, start, mid);
        node->right = buildBVH(indices, mid, end);
        return node;
    }

    void collideBVH(BVHNode *a, BVHNode *b)
    {
        if (!a || !b || !a->bounds.overlaps(b->bounds))
            return;

        if (a->isLeaf() && b->isLeaf())
        {
            if (a->bodyIndex != b->bodyIndex)
            {
                resolve(a->bodyIndex, b->bodyIndex);
            }
            return;
        }

        if (a->isLeaf())
        {
            collideBVH(a, b->left);
            collideBVH(a, b->right);
        }
        else if (b->isLeaf())
        {
            collideBVH(a->left, b);
            collideBVH(a->right, b);
        }
        else
        {
            collideBVH(a->left, b->left);
            collideBVH(a->left, b->right);
            collideBVH(a->right, b->left);
            collideBVH(a->right, b->right);
        }
    }

    void resolve(size_t i, size_t j)
    {
        const Body &b1 = bodies[i];
        const Body &b2 = bodies[j];

        Vec2 d = b2.pos - b1.pos;
        float r = b1.radius + b2.radius;

        if (d.mag_sq() > r * r)
            return;

        Vec2 v = b2.vel - b1.vel;
        float d_dot_v = d.dot(v);

        float m1 = b1.mass;
        float m2 = b2.mass;

        float weight1 = m2 / (m1 + m2);
        float weight2 = m1 / (m1 + m2);

        if (d_dot_v >= 0.0f && !(d == Vec2::zero()))
        {
            Vec2 tmp = d * (r / d.mag() - 1.0f);
            bodies[i].pos -= tmp * weight1;
            bodies[j].pos += tmp * weight2;
            return;
        }

        float v_sq = v.mag_sq();
        float d_sq = d.mag_sq();
        float r_sq = r * r;

        float discriminant = d_dot_v * d_dot_v - v_sq * (d_sq - r_sq);
        if (discriminant < 0.0f)
            discriminant = 0.0f;

        float t = (d_dot_v + std::sqrt(discriminant)) / v_sq;

        bodies[i].pos -= b1.vel * t;
        bodies[j].pos -= b2.vel * t;

        Vec2 new_d = bodies[j].pos - bodies[i].pos;
        float new_d_dot_v = new_d.dot(v);
        float new_d_sq = new_d.mag_sq();

        Vec2 tmp = new_d * (1.5f * new_d_dot_v / new_d_sq);
        Vec2 v1 = b1.vel + tmp * weight1;
        Vec2 v2 = b2.vel - tmp * weight2;

        bodies[i].vel = v1;
        bodies[j].vel = v2;
        bodies[i].pos += v1 * t;
        bodies[j].pos += v2 * t;
    }

    std::vector<Body> uniform_disc(size_t n)
    {
        std::mt19937 rng(0);
        std::uniform_real_distribution<float> dist(0.0f, 1.0f);

        float inner_radius = 100.0f;
        float outer_radius = std::sqrt(static_cast<float>(n)) * 300.7f;

        std::vector<Body> bodies;
        bodies.reserve(n);

        float m = 1e9f; // Central mass (e.g., black hole)
        bodies.push_back(Body(Vec2::zero(), Vec2::zero(), m, inner_radius));

        // const float TAU = 6.28318530718f; // Tau constant
        // const float TAU = (std::atan(1)*4)*2.0f;
        constexpr double TAU = std::numbers::pi * 2;

        // mass ranges and their probabilities here
        struct MassRange
        {
            float min_mass;
            float max_mass;
            float probability; 
        };

        std::vector<MassRange> massRanges = {
            {0.00005f, 0.8f, 0.825f}, // ~80–85% chance for masses between 0.005 and 1.0
            {1.2f, 2.5f, 0.125f},     // ~10–15% chance for masses between 1.0 and 3.0
            {5.0f, 50.0f, 0.025f}     // ~1–5% chance for masses between 3.0 and 10.0
        };

        // normalizes all of the probabilities to ensure they sum to 1.0
        float totalProbability = 0.0f;
        for (const auto &range : massRanges)
        {
            totalProbability += range.probability;
        }
        for (auto &range : massRanges)
        {
            range.probability /= totalProbability;
        }

        // and computes the cumulative probabilities for the range selection
        std::vector<float> cumulativeProbs;
        float cumulative = 0.0f;
        for (const auto &range : massRanges)
        {
            cumulative += range.probability;
            cumulativeProbs.push_back(cumulative);
        }

        while (bodies.size() < n)
        {
            double a = static_cast<double>(dist(rng)) * TAU;
            // normal
            // float sin_a = std::sin(a);
            // float cos_a = std::cos(a);


            /*chaos and cool stuff (m and outer radius will need to be tweaked)*/
            // double k = 6.0; // adjusting k changes the number of loops
            // float sin_a = std::sin(k * a);
            // float cos_a = std::sin(k * a) * std::cos(a);

            // lissajous curve
            // float sin_a = std::sin(a) * std::cos(a); // infinity
            // float cos_a = std::cos(a) * std::tan(a);

            // float denom = static_cast<double>(1.0f) + std::sin(a) * std::sin(a);
            // float sin_a = static_cast<double>(std::sqrt(2.0f)) * std::sin(a) * std::cos(a) / static_cast<double>(denom);
            // float cos_a = static_cast<double>(std::sqrt(2.0f)) * std::cos(a) / static_cast<double>(denom);

            double k = 7.0; // Number of petals
            // float o = outer_radius * sqrt(dist(rng));
            // float sin_a  = static_cast<double>(o) * std::cos(k * a) * std::cos(a);
            // float cos_a = static_cast<double>(o) * std::cos(k * a) * std::sin(a);

            // float sin_a = std::sin(a) / std::cos(a); // the grid / spiral
            // float cos_a = std::cos(a) / std::tan(a);

            // flip flop
            // float sin_a = std::sin(a) * std::cos(a * static_cast<double>(2.0f)); // double frequency spiral
            // float cos_a = std::cos(a) * std::sin(a * static_cast<double>(2.0f));

            // // // rose Lissajous curve
            // float sin_a = std::sin(a * 5.0) * std::cos(a); // 5 petal flower
            // float cos_a = std::cos(a * 5.0) * std::sin(a);

            // // box...
            // float sin_a = std::tanh(std::sin(a * 2.0));
            // float cos_a = std::tanh(std::cos(a * 2.0));

            // // star burst
            // float sin_a = std::sin(k * a) * std::pow(std::cos(a * 3.0), k);
            // float cos_a = std::cos(k * a) * std::pow(std::sin(a * 3.0), k);
            //using k
            float sin_a = std::sin(k * a) * std::pow(std::cos(a * 3.0), k);
            float cos_a = std::cos(k * a) * std::pow(std::sin(a * 3.0), k);

            // golden
            // float sin_a = std::sin(a) * std::pow(std::cosh(a * 3.0), 2);
            // float cos_a = std::cos(a) * std::pow(std::sinh(a * 3.0), 2);

            // // ?
            // float phase = std::fmod(frame * 0.01f, TAU);
            // float sin_a = std::tan(a + static_cast<double>(phase)) / std::cos(a / 0.5 + static_cast<double>(phase));
            // float cos_a = std::tanh(a + static_cast<double>(phase)) / std::sin(a / 0.5 + static_cast<double>(phase));

            // float sin_a = std::cos(a); //cool
            // float cos_a = std::tan(a);

            // float sin_a = std::tan(a); // interesting
            // float cos_a = std::cosh(a);

            // float sin_a = std::cosh(a); // interesting
            // float cos_a = std::tan(a);

            float t = inner_radius / outer_radius;
            float r = dist(rng) * (1.0f - t * t) + t * t;
            Vec2 pos(cos_a, sin_a);
            pos *= outer_radius * std::sqrt(r);
            Vec2 vel(sin_a, -cos_a);

            // will select a mass range based the abobe probabilities
            float randProb = dist(rng);
            size_t selectedRangeIndex = 0;
            for (size_t i = 0; i < cumulativeProbs.size(); ++i)
            {
                if (randProb <= cumulativeProbs[i])
                {
                    selectedRangeIndex = i;
                    break;
                }
            }
            const MassRange &selectedRange = massRanges[selectedRangeIndex];

            // assigns the mass within the selected range
            float mass = dist(rng) * (selectedRange.max_mass - selectedRange.min_mass) + selectedRange.min_mass;

            float radius = std::cbrt(mass); // and adjust radius based on mass

            bodies.push_back(Body(pos, vel, mass, radius));
        }

        // sort all of bodies by distance from the center
        std::sort(bodies.begin(), bodies.end(),
                  [](const Body &a, const Body &b)
                  {
                      return a.pos.mag_sq() < b.pos.mag_sq();
                  });

        // adjust velocities for circular orbits
        float total_mass = 0.0f;
        for (size_t i = 0; i < bodies.size(); ++i)
        {
            total_mass += bodies[i].mass;
            if (bodies[i].pos == Vec2::zero())
                continue;
            float v = std::sqrt(total_mass / bodies[i].pos.mag());
            bodies[i].vel *= v;
        }

        return bodies;
    }
};