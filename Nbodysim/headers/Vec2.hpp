#pragma once
#include <cmath>
#include <raylib.h>
#include <algorithm>
#include <thread>
#include <atomic>
#include <mutex>
#include <future>
#include <vector>

#ifdef __ARM_NEON
#include <arm_neon.h>
#elif defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif

struct alignas(16) Vec2
{

    float x, y;

    // constexpr Vec2() noexcept : x(0.0f), y(0.0f) {} //494ms of performance loss
    constexpr Vec2() noexcept = default;
    // constexpr Vec2(float x, float y) noexcept : x(x), y(y) {}
    constexpr Vec2(float x_, float y_) noexcept : x(x_), y(y_) {}
    inline constexpr Vec2(const Vec2 &other) noexcept : x(other.x), y(other.y) {}

    constexpr Vec2(Vec2 &&other) noexcept : x(other.x), y(other.y) {}

    static constexpr Vec2 zero() noexcept { return Vec2(0.0f, 0.0f); }
    static constexpr Vec2 one() noexcept { return Vec2(1.0f, 1.0f); }
    static constexpr Vec2 unit_x() noexcept { return Vec2(1.0f, 0.0f); }
    static constexpr Vec2 unit_y() noexcept { return Vec2(0.0f, 1.0f); }
    static constexpr Vec2 broadcast(float val) noexcept { return Vec2(val, val); }

#if defined(__ARM_NEON)
    Vec2 operator+(const Vec2 &other) const noexcept
    {
        Vec2 result;
        float32x2_t v1 = {x, y};
        float32x2_t v2 = {other.x, other.y};
        float32x2_t sum = vadd_f32(v1, v2);
        vst1_f32(&result.x, sum);
        return result;
    }
   inline Vec2 operator-(const Vec2 &other) const noexcept
    {
        Vec2 result;
        float32x2_t v1 = {x, y};
        float32x2_t v2 = {other.x, other.y};
        float32x2_t diff = vsub_f32(v1, v2);
        vst1_f32(&result.x, diff);
        return result;
    }

    inline Vec2 operator*(float scalar) const noexcept
    {
        float32x2_t a = vld1_f32(&x);
        float32x2_t s = vdup_n_f32(scalar); 
        float32x2_t result = vmul_f32(a, s);
        Vec2 ret;                           
        vst1_f32(&ret.x, result);
        return ret;
    }

    Vec2 operator/(float scalar) const noexcept
    {
        float32x2_t a = vld1_f32(&x);
        float32x2_t s = vdup_n_f32(1.0f / scalar);
        float32x2_t result = vmul_f32(a, s);
        Vec2 ret;
        vst1_f32(&ret.x, result);
        return ret;
    }

#elif defined(__x86_64__) || defined(_M_X64)
    Vec2 operator+(const Vec2 &other) const noexcept
    {
        __m128 a = _mm_setr_ps(x, y, 0, 0);
        __m128 b = _mm_setr_ps(other.x, other.y, 0, 0);
        __m128 r = _mm_add_ps(a, b);
        float res[4];
        _mm_store_ps(res, r);
        return Vec2(res[0], res[1]);
    }

    Vec2 operator-(const Vec2 &other) const noexcept
    {
        __m128 a = _mm_setr_ps(x, y, 0, 0);
        __m128 b = _mm_setr_ps(other.x, other.y, 0, 0);
        __m128 r = _mm_sub_ps(a, b);
        float res[4];
        _mm_store_ps(res, r);
        return Vec2(res[0], res[1]);
    }

    Vec2 operator*(float scalar) const noexcept
    {
        __m128 a = _mm_setr_ps(x, y, 0, 0);
        __m128 s = _mm_set1_ps(scalar);
        __m128 r = _mm_mul_ps(a, s);
        float res[4];
        _mm_store_ps(res, r);
        return Vec2(res[0], res[1]);
    }

    Vec2 operator/(float scalar) const noexcept
    {
        __m128 a = _mm_setr_ps(x, y, 0, 0);
        __m128 s = _mm_set1_ps(1.0f / scalar);
        __m128 r = _mm_mul_ps(a, s);
        float res[4];
        _mm_store_ps(res, r);
        return Vec2(res[0], res[1]);
    }
#else
    Vec2 operator+(const Vec2 &other) const noexcept
    {
        return Vec2(x + other.x, y + other.y);
    }
    Vec2 operator-(const Vec2 &other) const noexcept
    {
        return Vec2(x - other.x, y - other.y);
    }
    Vec2 operator*(float scalar) const noexcept
    {
        return Vec2(x * scalar, y * scalar);
    }
    Vec2 operator/(float scalar) const noexcept
    {
        return Vec2(x / scalar, y / scalar);
    }
#endif
    // std::span<const float> as_span() const noexcept
    // {
    //     return std::span<const float>(reinterpret_cast<const float *>(this), 2);
    // }

    inline Vec2 &operator+=(const Vec2 &other) noexcept
    {
        x += other.x; // a direct modification instead of *this = *this + other
        y += other.y;
        return *this;
    }

    Vec2 &operator-=(const Vec2 &other) noexcept
    {
        x -= other.x;
        y -= other.y;
        return *this;
    }
    Vec2 &operator*=(float scalar) noexcept
    {
        x *= scalar;
        y *= scalar;
        return *this;
    }

    Vec2 &operator/=(float scalar) noexcept
    {
        float inv = 1.0f / scalar;
        x *= inv;
        y *= inv;
        return *this;
    }

    constexpr bool operator==(const Vec2 &other) const noexcept
    {
        return x == other.x && y == other.y;
    }

    Vec2 &operator=(const Vec2 &other) noexcept
    {
        x = other.x;
        y = other.y; 
        return *this; 
    }

    Vec2 &operator=(Vec2 &&other) noexcept
    {
        x = other.x;
        y = other.y;
        return *this;
    }
    Vec2 mul_add(const Vec2 &mul, const Vec2 &add) const noexcept
    {
#ifdef __ARM_NEON
        float32x2_t v1 = {x, y};
        float32x2_t v2 = {mul.x, mul.y};
        float32x2_t v3 = {add.x, add.y};
        float32x2_t result = vfma_f32(v3, v1, v2);
        Vec2 ret;
        vst1_f32(&ret.x, result);
        return ret;
#else
        return Vec2(
            std::fma(x, mul.x, add.x),
            std::fma(y, mul.y, add.y));
#endif
    }

    float dot(const Vec2 &other) const noexcept
    {
#ifdef __x86_64__
        __m128 a = _mm_setr_ps(x, y, 0, 0);
        __m128 b = _mm_setr_ps(other.x, other.y, 0, 0);
        __m128 r = _mm_dp_ps(a, b, 0x31);
        float res;
        _mm_store_ss(&res, r);
        return res;
#else
        return x * other.x + y * other.y;
#endif
    }

    float mag_sq() const noexcept
    {
        return x * x + y * y;
    }

    float mag() const noexcept
    {
        return std::sqrt(mag_sq());
    }

    Vec2 &normalize() noexcept
    {
        float m = mag();
        if (m > 0.0f)
        {
            x /= m;
            x /= m;
            y /= m;
        }
        return *this;
    }

    [[nodiscard]] Vec2 normalized() const noexcept
    {
        Vec2 result = *this;
        float m = mag();
        if (m > 0.0f)
        {
            result.x /= m;
            result.y /= m;
        }
        return result;
    }

    Vec2 reflected(const Vec2 &normal) const noexcept
    {
        return *this - normal * (2.0f * dot(normal));
    }

    Vec2 clamped(const Vec2 &min, const Vec2 &max) const noexcept
    {
        return Vec2(
            std::clamp(x, min.x, max.x),
            std::clamp(y, min.y, max.y));
    }

    float component_max() const noexcept
    {
        return std::max(x, y);
    }

    float component_min() const noexcept
    {
        return std::min(x, y);
    }

    Vec2 abs() const noexcept
    {
        return Vec2(std::abs(x), std::abs(y));
    }

    Vec2 min_by_component(const Vec2 &other) const noexcept
    {
        return Vec2(std::min(x, other.x), std::min(y, other.y));
    }

    Vec2 max_by_component(const Vec2 &other) const noexcept
    {
        return Vec2(std::max(x, other.x), std::max(y, other.y));
    }

    float &operator[](size_t index) noexcept
    {
        return (index == 0) ? x : y;
    }

    const float &operator[](size_t index) const noexcept
    {
        return (index == 0) ? x : y;
    }

    static Vec2 min(const Vec2 &a, const Vec2 &b) noexcept
    {
        return Vec2(std::min(a.x, b.x), std::min(a.y, b.y));
    }

    static Vec2 max(const Vec2 &a, const Vec2 &b) noexcept
    {
        return Vec2(std::max(a.x, b.x), std::max(a.y, b.y));
    }

    operator Vector2() const noexcept { return Vector2{x, y}; }
    static Vec2 from_vector2(Vector2 v) noexcept { return Vec2(v.x, v.y); }
    Vector2 toVector2() const { return Vector2{x, y}; }

    static void batch_accumulate(Vec2 *acc, const Vec2 *pos,
                                 const Vec2 *com, const float *mass,
                                 size_t n) noexcept
    {
#ifdef __ARM_NEON
        for (size_t i = 0; i < n; i += 2)
        {
            float32x2_t p = vld1_f32(&pos[i].x);
            float32x2_t c = vld1_f32(&com[i].x);
            float32x2_t m = vld1_f32(&mass[i]);

            float32x2_t d = vsub_f32(c, p);
            float32x2_t d_sq = vmul_f32(d, d);

            uint32x2_t mask = vcgt_f32(d_sq, vdup_n_f32(1e-8f));
            d_sq = vbsl_f32(mask, d_sq, vdup_n_f32(1e-8f));

            float32x2_t inv_d = vrsqrte_f32(d_sq);
            inv_d = vmul_f32(inv_d, vrsqrts_f32(vmul_f32(d_sq, inv_d), inv_d));

            float32x2_t force = vmul_f32(m, vmul_f32(inv_d, vmul_f32(inv_d, inv_d)));
            float32x2_t result = vmul_f32(d, force);

            float32x2_t a = vld1_f32(&acc[i].x);
            a = vadd_f32(a, result);
            vst1_f32(&acc[i].x, a);
        }
#else
        for (size_t i = 0; i < n; ++i)
        {
            Vec2 d = com[i] - pos[i];
            float d_sq = d.mag_sq();
            if (d_sq > 0)
            {
                float inv_dist = 1.0f / std::sqrt(d_sq);
                acc[i] += d * (mass[i] * inv_dist * inv_dist * inv_dist);
            }
        }
#endif
    }

};

inline Vec2 operator*(float scalar, const Vec2 &vec) noexcept
{
    return vec * scalar;
}