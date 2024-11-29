#pragma once
#include "Body.hpp"
#include "Quadtree.hpp"
#include "AABB.hpp"
#include <vector>
#include <random>
#include <numeric>
#include <cmath>
#include <numbers>
#include <unordered_map>
#include <atomic>
#include <thread>
#include <future>
#include <execution>

extern std::atomic<float> SIMULATION_DT;


struct SpatialGrid
{
    static constexpr float CELL_SIZE = 600.0f;
    static constexpr size_t ESTIMATED_BODIES_PER_CELL = 16;

    struct Cell
    {
        std::vector<size_t> bodies;
        Cell() { bodies.reserve(ESTIMATED_BODIES_PER_CELL); }
    };

    std::unordered_map<size_t, Cell> cells;

    static size_t hash_position(int x, int y)
    {
        return ((x * 92837111) ^ (y * 689287499)) * 15485863;
    }
};


struct SweepEntry
{
    float value;    // min_x or max_x of a body's AABB
    size_t bodyIdx; // Index of the body in the bodies vector
    bool isEnd;     // true if it's the max_x value (end of the interval) then false if min_x

    bool operator<(const SweepEntry &other) const
    {
        return value < other.value || (value == other.value && isEnd < other.isEnd);
    }
};

class Simulation
{
public:
    float dt;
    size_t frame;
    std::vector<Body> bodies;
    Quadtree quadtree;
    std::vector<size_t> bodyIndices;

    Simulation()
        : frame(0), quadtree(1.0f, 1.0f, 16)
    {
        size_t n = 100000;
        bodies.reserve(n);
        bodyIndices.reserve(n);
        bodies = uniform_disc(n);
    }

    void step()
    {
        float current_dt = SIMULATION_DT.load();
        iterate(current_dt);
        //  Body::batch_update_parallel(bodies, current_dt);
        collide();
        attract();
        ++frame;

        
    }

private:
    void iterate(float dt)
    {

        for (auto &body : bodies)
        {
            body.update(dt);
        }
    }

   
    // significant bottleneck
    //  void attract()
    //  {
    //      quadtree.build(bodies);

    //     for (size_t i = 0; i < bodies.size(); i++)
    //     {
    //         bodies[i].acc = quadtree.acc(bodies[i].pos, bodies);
    //     }
    // }

    void attract()
    {
        quadtree.build(bodies);

        const unsigned int num_threads = std::thread::hardware_concurrency() ? std::thread::hardware_concurrency() : 4;
        std::vector<std::future<void>> futures;

        size_t chunk_size = bodies.size() / num_threads;

        for (unsigned int t = 0; t < num_threads; ++t)
        {
            size_t start = t * chunk_size;
            size_t end = (t == num_threads - 1) ? bodies.size() : start + chunk_size;

            futures.emplace_back(std::async(std::launch::async, [&, start, end]()
                                            {
            for (size_t i = start; i < end; ++i) {
                bodies[i].acc = quadtree.acc(bodies[i].pos, bodies); 
            } }));
        }

        for (auto &fut : futures)
        {
            fut.get();
        }
    }

    void collide()
    {
        if (bodies.empty())
            return;

        SpatialGrid grid;
        std::vector<std::pair<size_t, size_t>> broadPhasePairs;

        for (size_t i = 0; i < bodies.size(); i++)
        {
            const Body &body = bodies[i];
            Vec2 r(body.radius, body.radius);
            AABB bounds = {body.pos - r, body.pos + r};

            int minX = static_cast<int>(bounds.min.x / SpatialGrid::CELL_SIZE);
            int maxX = static_cast<int>(bounds.max.x / SpatialGrid::CELL_SIZE);
            int minY = static_cast<int>(bounds.min.y / SpatialGrid::CELL_SIZE);
            int maxY = static_cast<int>(bounds.max.y / SpatialGrid::CELL_SIZE);

            for (int y = minY; y <= maxY; y++)
            {
                for (int x = minX; x <= maxX; x++)
                {
                    size_t hash = SpatialGrid::hash_position(x, y);
                    grid.cells[hash].bodies.push_back(i);
                }
            }
        }

        for (const auto &[hash, cell] : grid.cells)
        {
            if (cell.bodies.size() > 1)
            {
                std::vector<SweepEntry> sweepList;
                sweepList.reserve(cell.bodies.size() * 2);

                for (size_t idx : cell.bodies)
                {
                    const Body &body = bodies[idx];
                    Vec2 r(body.radius, body.radius);
                    float minX = body.pos.x - r.x;
                    float maxX = body.pos.x + r.x;

                    sweepList.push_back({minX, idx, false});
                    sweepList.push_back({maxX, idx, true});
                }

                std::sort(sweepList.begin(), sweepList.end());

                std::vector<size_t> active;

                for (const auto &entry : sweepList)
                {
                    if (!entry.isEnd)
                    {
                        for (size_t activeIdx : active)
                        {
                            broadPhasePairs.emplace_back(activeIdx, entry.bodyIdx);
                        }
                        active.push_back(entry.bodyIdx);
                    }
                    else
                    {
                        active.erase(std::remove(active.begin(), active.end(),
                                                 entry.bodyIdx),
                                     active.end());
                    }
                }
            }
        }


        for (const auto &[i, j] : broadPhasePairs)
        {
            resolve(i, j);
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

        float inner_radius = 200.0f;
        float outer_radius = std::sqrt(static_cast<float>(n)) * 300.7f;

        std::vector<Body> bodies;
        bodies.reserve(n);

        float m = 1e9f; // Central mass (e.g., black hole)
        bodies.push_back(Body(Vec2::zero(), Vec2::zero(), m, inner_radius));

        // const float TAU = 6.28318530718f; // Tau constant
        // const float TAU = (std::atan(1)*4)*2.0f;
        // constexpr double TAU = std::numbers::pi * 4;

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
        // lorenz attractor parameters
        const float sigma = 10.0f;
        const float rho = 28.0f;
        const float beta = 8.0f / 3.0f;

        float x = 0.1f;
        float y = 0.0f;
        float z = 0.0f;

        while (bodies.size() < n)
        {
            // double a = static_cast<double>(dist(rng)) * TAU;
            // double k = 7.0; // adjusting k changes the number of loops and other fun things

            // normal
            // float sin_a = std::sin(a);
            // float cos_a = std::cos(a);

            /*chaos and cool stuff (m and outer radius will need to be tweaked)*/
            // float sin_a = std::sin(k * a);
            // float cos_a = std::sin(k * a) * std::cos(a);

            // lissajous curve
            // float sin_a = std::sin(a) * std::cos(a); // infinity
            // float cos_a = std::cos(a) * std::tan(a);

            // float denom = static_cast<double>(1.0f) + std::sin(a) * std::sin(a);
            // float sin_a = static_cast<double>(std::sqrt(2.0f)) * std::sin(a) * std::cos(a) / static_cast<double>(denom);
            // float cos_a = static_cast<double>(std::sqrt(2.0f)) * std::cos(a) / static_cast<double>(denom);

            // double k = 7.0; // Number of petals
            // float o = outer_radius * sqrt(dist(rng));
            // float sin_a  = static_cast<double>(o) * std::cos(k * a) * std::cos(a);
            // float cos_a = static_cast<double>(o) * std::cos(k * a) * std::sin(a);

            // float sin_a = std::sin(a * k) / std::cos(a); // the grid / spiral
            // float cos_a = std::cos(a * k) / std::tan(a);

            // Heart shape parametric equations <3
            // float f = static_cast<float>(a) * 2.0f; // Force float precision
            // float scale = 1.2f;

            // // Break down calculations for clarity and precision
            // float sin_cube = std::powf(std::sinf(f), 3.0f);
            // float sin_term = 16.0f * sin_cube;

            // // Cosine terms with decreasing amplitudes
            // float cos_term1 = 13.0f * std::cosf(f);
            // float cos_term2 = 5.0f * std::cosf(2.0f * f);
            // float cos_term3 = 2.0f * std::cosf(3.0f * f);
            // float cos_term4 = std::cosf(4.0f * f);

            // // Final parametric equations
            // float sin_a = scale * sin_term;
            // float cos_a = scale * (cos_term1 - cos_term2 - cos_term3 - cos_term4);

            // Heart shape parametric equations - vertical orientation
            // float f = static_cast<float>(a) * 2.0f;
            // float scale = 1.2f;

            // // Use powf and single-precision trig functions consistently
            // float sin_v = std::sinf(f);
            // float sin_cube = sin_v * sin_v * sin_v;  // Instead of powf
            // float sin_term = 16.0f * sin_cube;

            // // Cosine terms using single-precision
            // float cos_term1 = 13.0f * std::cosf(f);
            // float cos_term2 = 5.0f * std::cosf(2.0f * f);
            // float cos_term3 = 2.0f * std::cosf(3.0f * f);
            // float cos_term4 = std::cosf(4.0f * f);

            // // Final parametric equations
            // float cos_a = scale * sin_term;  // x coordinate
            // float sin_a = -scale * (cos_term1 - cos_term2 - cos_term3 - cos_term4);  // y coordinate
            // flip flop
            // float sin_a = std::sin(k*a) * std::cos(a * static_cast<double>(2.0f)); // double frequency spiral
            // float cos_a = std::cos(k*a) * std::sin(a * static_cast<double>(2.0f));

            // // // rose Lissajous curve
            // float sin_a = std::sin(a * 5.0) * std::cos(a); // 5 petal flower
            // float cos_a = std::cos(a * 5.0) * std::sin(a);

            // // box...
            // float sin_a = std::tanh(std::sin(k* a * 2.0));
            // float cos_a = std::tanh(std::cos(k* a * 2.0));

            // // star burst
            // float sin_a = std::sin(k * a) * std::pow(std::cos(a * 3.0), k);
            // float cos_a = std::cos(k * a) * std::pow(std::sin(a * 3.0), k);
            // using k
            // float sin_a = std::sin(k * a) * std::pow(std::cos(a * 3.0), k);
            // float cos_a = std::cos(k * a) * std::pow(std::sin(a * 3.0), k);

            // The mural
            // Create oscillating patterns using frame count
            // float phase = std::fmod(frame * 0.01f, TAU);  // Smooth oscillation over time
            // float blend = (std::asinh(phase) + 1.0f) * 0.5f;  // Oscillates between 0 and 1

            // float sin_a = blend * (std::tanf(k * a) * std::powf(std::cosf(a * 3.0), k)) +
            //               (1.0f + blend) * (std::tanhf(k * a) * std::powf(std::sinf(a * 3.0), k));

            // float cos_a = blend * (std::tanf(k * a) * std::powf(std::sinf(a * 3.0), k)) +
            //               (1.0f + blend) * (std::tanhf(k * a) * std::powf(std::cosf(a * 3.0), k));

            // golden
            // float sin_a = std::sin(a) * std::pow(std::cosh(a * 3.0), 2);
            // float cos_a = std::cos(a) * std::pow(std::sinh(a * 3.0), 2);

            // // ?
            // float phase = std::fmod(frame * 0.01f, TAU);
            // float sin_a = std::tan(a + static_cast<double>(phase)) / std::cos(k*a / 0.5 + static_cast<double>(phase));
            // float cos_a = std::tanh(a + static_cast<double>(phase)) / std::sin(k*a / 0.5 + static_cast<double>(phase));

            // float sin_a = std::cos(a); //cool
            // float cos_a = std::tan(a);

            // double b = static_cast<double>(dist(rng)) * TAU; // Assuming y is also a random variable

            // // // Equation sin(x^2 + y^2) = cos(x)
            // float sin_a = std::sin(a * a + b * b); // sin(x^2 + y^2)
            // float cos_a = std::cos(a);             // cos(x)

            // // Combine the results
            // float result = sin_a - cos_a; // T

            const float dt = 0.01f; // Integration timestep for Lorenz

            // Lorenz equations
            float dx = sigma * (y - x);
            float dy = x * (rho - z) - y;
            float dz = x * y - beta * z;

            x += dx * dt;
            y += dy * dt;
            z += dz * dt;

            float scale = outer_radius / 10.0f; 
            Vec2 pos(x * scale, y * scale);

            Vec2 vel(-pos.y, pos.x);
            vel.normalize();

            // lemniscate of Bernoulli
            // float sin_a = std::cos(a) / (1.0 + std::sin(a) * std::sin(a));
            // float cos_a = std::cos(a) * std::sin(a) / (1.0 + std::sin(a) * std::sin(a));

            // fermat's spiral
            // float c = 0.5f; // controls spacing
            // float r1 = c * std::sqrtf(a);
            // float sin_a = r1 * std::sinf(a);
            // float cos_a = r1 * std::cosf(a);
            // hypocycloid (4-cusped) (eye)
            // float k = 4.0f;
            // float sin_a = (k - 1.0f) * std::cosf(static_cast<float>(a)) + std::cosf((k - 1.0f) * static_cast<float>(a));
            // float cos_a = (k - 1.0f) * std::sinf(static_cast<float>(a)) - std::sinf((k - 1.0f) * static_cast<float>(a));
            //
            // float k = 5.0f;
            // float sin_a = (k + 1.0f) * std::cosf(static_cast<float>(a)) + std::cosf((k + 1.0f) * static_cast<float>(a));
            // float cos_a = (k + 1.0f) * std::sinf(static_cast<float>(a)) - std::sinf((k + 1.0f) * static_cast<float>(a));

            // float t = inner_radius / outer_radius;
            // float r = dist(rng) * (1.0f - t * t) + t * t;
            // Vec2 pos(cos_a, sin_a); // Use x and y for position
            // pos *= outer_radius * std::sqrt(r);
            // Vec2 vel(sin_a, -cos_a);

            // will select a mass range based the above probabilities
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

            float mass = dist(rng) * (selectedRange.max_mass - selectedRange.min_mass) + selectedRange.min_mass;

            float radius = std::cbrt(mass); 

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