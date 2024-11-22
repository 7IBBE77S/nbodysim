#pragma once
#include <cmath>
#include <raylib.h>

#ifdef __ARM_NEON
#include <arm_neon.h>
#elif defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif
// #include <raymath.h>

struct alignas(16) Vec2
{
    float x, y;

    Vec2() : x(0.0f), y(0.0f) {}
    Vec2(float x, float y) : x(x), y(y) {}

    static Vec2 zero() { return Vec2(0.0f, 0.0f); }
    static Vec2 one() { return Vec2(1.0f, 1.0f); }

    Vec2 operator+(const Vec2 &other) const { return Vec2(x + other.x, y + other.y); }
    Vec2 operator-(const Vec2 &other) const { return Vec2(x - other.x, y - other.y); }
    Vec2 operator*(float scalar) const { return Vec2(x * scalar, y * scalar); }
    Vec2 operator/(float scalar) const { return Vec2(x / scalar, y / scalar); }

    Vec2 &operator+=(const Vec2 &other)
    {
        x += other.x;
        y += other.y;
        return *this;
    }

    Vec2 &operator-=(const Vec2 &other)
    {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    Vec2 &operator*=(float scalar)
    {
        x *= scalar;
        y *= scalar;
        return *this;
    }

    Vec2 &operator/=(float scalar)
    {
        x /= scalar;
        y /= scalar;
        return *this;
    }

    bool operator==(const Vec2 &other) const
    {
        return x == other.x && y == other.y;
    }

    float dot(const Vec2 &other) const
    {
        return x * other.x + y * other.y;
    }

    float mag_sq() const
    {
        return x * x + y * y;
    }

    float mag() const
    {
        return std::sqrt(mag_sq());
    }

    // Converts to the raylib Vector2
    // operator Vector2() const {
    //     return { x, y };
    // }
    operator Vector2() const { return Vector2{x, y}; }

    void normalize()
    {
        float mag = std::sqrt(x * x + y * y);
        if (mag > 0.0f)
        {
            x /= mag;
            y /= mag;
        }
    }

    //UNUSED
 #ifdef __ARM_NEON
    // ARM NEON version
    static void batch_accumulate(
        Vec2 *acc,
        const Vec2 *pos,
        const Vec2 *center_of_mass,
        const float *mass,
        size_t n)
    {
        for (size_t i = 0; i < n; i += 4)
        {
            float32x4_t px = vld1q_f32(&pos[i].x);
            float32x4_t py = vld1q_f32(&pos[i].y);
            float32x4_t cx = vld1q_f32(&center_of_mass[i].x);
            float32x4_t cy = vld1q_f32(&center_of_mass[i].y);
            float32x4_t m = vld1q_f32(&mass[i]);

            // calculate the displacement
            float32x4_t dx = vsubq_f32(cx, px);
            float32x4_t dy = vsubq_f32(cy, py);

            // calculate distance squared
            float32x4_t d_sq = vaddq_f32(
                vmulq_f32(dx, dx),
                vmulq_f32(dy, dy));

            // calcs force
            float32x4_t force = vdivq_f32(m,
                                          vmulq_f32(d_sq, vsqrtq_f32(d_sq)));

            // accumulate acceleration
            vst1q_f32(&acc[i].x, vmulq_f32(dx, force));
            vst1q_f32(&acc[i].y, vmulq_f32(dy, force));
        }
    }
#else
    // fallback non SIMD version
    static void batch_accumulate(
        Vec2 *acc,
        const Vec2 *pos,
        const Vec2 *center_of_mass,
        const float *mass,
        size_t n)
    {
        for (size_t i = 0; i < n; i++)
        {
            Vec2 d = center_of_mass[i] - pos[i];
            float d_sq = d.mag_sq();
            float force = mass[i] / (d_sq * std::sqrt(d_sq));
            acc[i] += d * force;
        }
    }
#endif
};
