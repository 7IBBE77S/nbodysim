#pragma once
#include "Vec2.hpp"
#include <vector>
#include <future>
#include <thread>
struct alignas(16) Body
{
    Vec2 pos;     // 8 bytes
    Vec2 vel;     // 8 bytes
    Vec2 acc;     // 8 bytes
    float mass;   // 4 bytes
    float radius; // 4 bytes
                  // total: 32 bytes cache aligned

    Body() = default;
    Body(Vec2 pos, Vec2 vel, float mass, float radius)
        : pos(pos), vel(vel), acc(Vec2::zero()), mass(mass), radius(radius) {}

    static Body at_rest(Vec2 pos, float mass, float radius)
    {
        return Body(pos, Vec2::zero(), mass, radius);
    }

    static Body with_velocity(Vec2 pos, Vec2 vel, float mass, float radius)
    {
        return Body(pos, vel, mass, radius);
    }
    static unsigned int get_hardware_threads()
    {
        return std::thread::hardware_concurrency();
    }

    // physics update
   inline void update(float dt)
    {
        vel += acc * dt; // 152ms of performance loss
        pos += vel * dt; //230ms of performance loss
    }

 

    // Batch update test function
    static void batch_update_parallel(std::vector<Body> &bodies, float dt)
    {
        unsigned int num_threads = get_hardware_threads();
        std::vector<std::future<void>> futures;

        size_t chunk_size = bodies.size() / num_threads;

        for (unsigned int t = 0; t < num_threads; ++t)
        {
            size_t start = t * chunk_size;
            size_t end = (t == num_threads - 1) ? bodies.size() : (t + 1) * chunk_size;

            futures.push_back(std::async(std::launch::async, [&bodies, dt, start, end]()
                                         {
                for (size_t i = start; i < end; ++i) {
                    bodies[i].update(dt);
                } }));
        }

        for (auto &future : futures)
        {
            future.wait();
        }
    }

    bool collides_with(const Body &other) const
    {
        Vec2 d = other.pos - pos;
        float r = radius + other.radius;
        return d.mag_sq() <= r * r;
    }

    void resolve_collision(Body &other)
    {
        Vec2 d = other.pos - pos;
        float r = radius + other.radius;

        if (d.mag_sq() > r * r)
            return;

        Vec2 n = d.normalized();
        Vec2 v = other.vel - vel;

        float j = -(1.0f + 0.5f) * v.dot(n);
        j /= 1.0f / mass + 1.0f / other.mass;

        vel -= n * (j / mass);
        other.vel += n * (j / other.mass);

        float overlap = r - d.mag();
        Vec2 separation = n * (overlap * 0.5f);
        pos -= separation;
        other.pos += separation;
    }

    float kinetic_energy() const
    {
        return 0.5f * mass * vel.mag_sq();
    }

    Vec2 momentum() const
    {
        return vel * mass;
    }
};

// class BodySystem
// {
// private:
//     // Hot data (accessed every frame during force calculation)
//     struct HotData
//     {
//         std::vector<float> x; // Positions split for better SIMD
//         std::vector<float> y;
//         std::vector<float> mass;
//     };

//     // warm data (accessed during integration)
//     struct WarmData
//     {
//         std::vector<float> vx; // Velocities split
//         std::vector<float> vy;
//         std::vector<float> ax; // Accelerations split
//         std::vector<float> ay;
//     };

//     // Cold data (rarely accessed)
//     struct ColdData
//     {
//         std::vector<float> radius;
//     };

//     HotData hot;
//     WarmData warm;
//     ColdData cold;
//     size_t count = 0;

// public:
//     BodySystem() : count(0) {}

//     BodySystem(const BodySystem &other)
//     {
//         hot = other.hot;
//         warm = other.warm;
//         cold = other.cold;
//         count = other.count;
//     }

//     BodySystem(BodySystem &&other) noexcept
//     {
//         hot = std::move(other.hot);
//         warm = std::move(other.warm);
//         cold = std::move(other.cold);
//         count = other.count;
//         other.count = 0;
//     }
//     void reserve(size_t n)
//     {
//         hot.x.reserve(n);
//         hot.y.reserve(n);
//         hot.mass.reserve(n);
//         warm.vx.reserve(n);
//         warm.vy.reserve(n);
//         warm.ax.reserve(n);
//         warm.ay.reserve(n);
//         cold.radius.reserve(n);
//     }
//     
//     void swap_bodies(size_t i, size_t j)
//     {
//         std::swap(hot.x[i], hot.x[j]);
//         std::swap(hot.y[i], hot.y[j]);
//         std::swap(hot.mass[i], hot.mass[j]);
//         std::swap(warm.vx[i], warm.vx[j]);
//         std::swap(warm.vy[i], warm.vy[j]);
//         std::swap(warm.ax[i], warm.ax[j]);
//         std::swap(warm.ay[i], warm.ay[j]);
//         std::swap(cold.radius[i], cold.radius[j]);
//     }

//     void add_body(const Body &body)
//     {
//         hot.x.push_back(body.pos.x);
//         hot.y.push_back(body.pos.y);
//         hot.mass.push_back(body.mass);
//         warm.vx.push_back(body.vel.x);
//         warm.vy.push_back(body.vel.y);
//         warm.ax.push_back(body.acc.x);
//         warm.ay.push_back(body.acc.y);
//         cold.radius.push_back(body.radius);
//         count++;
//     }

//     void update(float dt)
//     {
//         if (count == 0)
//             return; // Early exit if empty

//         for (size_t i = 0; i < count; i++)
//         {
//             // Check array bounds
//             if (i >= hot.x.size() || i >= hot.y.size() ||
//                 i >= warm.vx.size() || i >= warm.vy.size())
//             {
//                 break;
//             }

//             warm.vx[i] += warm.ax[i] * dt;
//             warm.vy[i] += warm.ay[i] * dt;
//             hot.x[i] += warm.vx[i] * dt;
//             hot.y[i] += warm.vy[i] * dt;
//         }
//     }

//     // Getters/setters
//     Vec2 get_position(size_t i) const
//     {
//         return Vec2(hot.x[i], hot.y[i]);
//     }

//     void set_position(size_t i, const Vec2 &pos)
//     {
//         hot.x[i] = pos.x;
//         hot.y[i] = pos.y;
//     }

//     float get_mass(size_t i) const
//     {
//         return hot.mass[i];
//     }

//     void set_acceleration(size_t i, const Vec2 &acc)
//     {
//         warm.ax[i] = acc.x;
//         warm.ay[i] = acc.y;
//     }

//     Vec2 get_velocity(size_t i) const
//     {
//         return Vec2(warm.vx[i], warm.vy[i]);
//     }

//     float get_radius(size_t i) const
//     {
//         return cold.radius[i];
//     }

//     void set_velocity(size_t i, const Vec2 &vel)
//     {
//         warm.vx[i] = vel.x;
//         warm.vy[i] = vel.y;
//     }

//     size_t size() const { return count; }

//     BodySystem &operator=(BodySystem &&other) noexcept
//     {
//         if (this != &other)
//         {
//             hot = std::move(other.hot);
//             warm = std::move(other.warm);
//             cold = std::move(other.cold);
//             count = other.count;
//             other.count = 0;
//         }
//         return *this;
//     }
//     BodySystem& operator=(const BodySystem& other) {
//     if (this != &other) {
//         hot = other.hot;
//         warm = other.warm;
//         cold = other.cold;
//         count = other.count;
//     }
//     return *this;
// }
// };