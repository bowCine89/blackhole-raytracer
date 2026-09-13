// core.hpp -- small math types, RNG and sampling helpers.
#pragma once
#include <cmath>
#include <cstdint>
#include <algorithm>

using Real = double;   // geodesic integration wants the mantissa; keep it double.

constexpr Real PI     = 3.14159265358979323846;
constexpr Real TWO_PI = 6.28318530717958647693;
constexpr Real INV_PI  = 0.31830988618379067154;
constexpr Real HALF_PI = 1.57079632679489661923;

struct Vec3 {
    Real x = 0, y = 0, z = 0;
    Vec3() = default;
    Vec3(Real a, Real b, Real c) : x(a), y(b), z(c) {}
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(Real s)        const { return {x * s, y * s, z * s}; }
    Vec3 operator-()              const { return {-x, -y, -z}; }
    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
};
inline Real dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline Real length(const Vec3& a) { return std::sqrt(dot(a, a)); }
inline Vec3 normalize(const Vec3& a) { Real l = length(a); return l > 0 ? a * (1 / l) : a; }

// A 4-vector in a coordinate basis (index order t, r, theta, phi).
struct Vec4 {
    Real v[4] = {0, 0, 0, 0};
    Real& operator[](int i) { return v[i]; }
    Real  operator[](int i) const { return v[i]; }
};

// Diagonal-plus-(t,phi) metric of a stationary axisymmetric spacetime.
// Stored as a full 4x4 so the Gram-Schmidt routine stays generic.
struct Metric {
    Real g[4][4] = {};
    Real dot(const Vec4& a, const Vec4& b) const {
        Real s = 0;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                if (g[i][j] != 0) s += g[i][j] * a[i] * b[j];
        return s;
    }
    Vec4 lower(const Vec4& u) const {   // u^mu -> u_mu
        Vec4 d;
        for (int i = 0; i < 4; ++i) {
            Real s = 0;
            for (int j = 0; j < 4; ++j) s += g[i][j] * u[j];
            d[i] = s;
        }
        return d;
    }
};

// ---------------------------------------------------------------------------
// PCG32 -- tiny, fast, and statistically sound enough for Monte Carlo.
// ---------------------------------------------------------------------------
struct Rng {
    uint64_t state = 0x853c49e6748fea9bULL, inc = 0xda3e39cb94b95bdbULL;

    explicit Rng(uint64_t seq = 1, uint64_t seed = 0x4d595df4d0f33173ULL) {
        state = 0u; inc = (seq << 1u) | 1u;
        nextU32(); state += seed; nextU32();
    }
    uint32_t nextU32() {
        uint64_t old = state;
        state = old * 6364136223846793005ULL + inc;
        uint32_t xorshifted = uint32_t(((old >> 18u) ^ old) >> 27u);
        uint32_t rot = uint32_t(old >> 59u);
        return (xorshifted >> rot) | (xorshifted << ((~rot + 1u) & 31));
    }
    // Uniform in [0,1).
    Real uniform() { return (nextU32() >> 8) * (1.0 / 16777216.0); }
};

// Deterministic integer hash (used by the procedural star field).
inline uint32_t hashU32(uint32_t x) {
    x ^= x >> 16; x *= 0x7feb352dU;
    x ^= x >> 15; x *= 0x846ca68bU;
    x ^= x >> 16; return x;
}
inline uint32_t hashCombine(uint32_t a, uint32_t b, uint32_t c) {
    return hashU32(a ^ hashU32(b * 0x9e3779b9U ^ hashU32(c * 0x85ebca6bU)));
}
inline Real hashFloat(uint32_t h) { return (h >> 8) * (1.0 / 16777216.0); }

// Cosine-weighted hemisphere sample about +z (concentric disk mapping).
inline Vec3 cosineHemisphere(Real u1, Real u2) {
    Real r = std::sqrt(u1);
    Real phi = TWO_PI * u2;
    return {r * std::cos(phi), r * std::sin(phi), std::sqrt(std::max<Real>(0, 1 - u1))};
}

inline Real clampf(Real v, Real lo, Real hi) { return v < lo ? lo : (v > hi ? hi : v); }
