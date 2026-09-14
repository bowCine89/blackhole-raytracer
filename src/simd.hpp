// simd.hpp -- a thin 8-wide double vector layer, plus the transcendentals the
// geodesic right-hand side needs.
//
// Built on GCC/Clang vector extensions rather than raw intrinsics: the operators
// read like scalar code, and the same source compiles to one AVX-512 register
// per vector, two AVX2 registers, or four SSE2 ones, so nothing here requires
// AVX-512 to be correct -- only to be fast.
//
// The reason this file exists at all: a single ray's Dormand-Prince stages are
// strictly sequential, so scalar tracing is latency-bound at ~850 cycles/step
// with the FPU mostly idle.  Putting eight *independent rays* in eight lanes
// turns that latency into throughput.  The blocker was sin/cos -- there is no
// vector libm here -- which is what vsincos below supplies.
#pragma once
#include "core.hpp"
#include <cstdint>

#if !(defined(__clang__) || defined(__GNUC__))
#  error "simd.hpp needs GCC/Clang vector extensions"
#endif

// Lane count is a build-time knob because the right width is a measurement,
// not a given: see the ISA table in the README.
#ifndef KERR_LANES
#  define KERR_LANES 8
#endif
static constexpr int LANES = KERR_LANES;

using vd = double    __attribute__((vector_size(KERR_LANES * 8)));
using vi = long long __attribute__((vector_size(KERR_LANES * 8)));

static inline vd vsplat(double x) { vd r; for (int i = 0; i < LANES; ++i) r[i] = x; return r; }
static inline vi isplat(long long x) { vi r; for (int i = 0; i < LANES; ++i) r[i] = x; return r; }

static inline vd vbits(vi x) { return __builtin_bit_cast(vd, x); }
static inline vi ibits(vd x) { return __builtin_bit_cast(vi, x); }

// Lanewise select.  Comparison results are all-ones or all-zeros per lane, so
// this is a bit blend and stays correct for NaN and infinity operands.
static inline vd vsel(vi mask, vd a, vd b) {
    return vbits((mask & ibits(a)) | (~mask & ibits(b)));
}
static inline vi vsel(vi mask, vi a, vi b) { return (mask & a) | (~mask & b); }

static inline vd vabs(vd x)  { return vbits(ibits(x) & isplat(0x7fffffffffffffffLL)); }
static inline vd vneg(vd x)  { return vbits(ibits(x) ^ isplat(int64_t(0x8000000000000000ULL))); }
static inline vd vmin(vd a, vd b) { return vsel(a < b, a, b); }
static inline vd vmax(vd a, vd b) { return vsel(a > b, a, b); }
// __builtin_elementwise_sqrt is a Clang extension with no GCC counterpart, so
// GCC gets a lanewise loop instead.  Both spellings compile to one vsqrtpd per
// register at -O3 -- the loop is a different way of writing it, not a slow path.
#if defined(__clang__)
static inline vd vsqrt(vd x) { return __builtin_elementwise_sqrt(x); }
#else
static inline vd vsqrt(vd x) {
    vd r;
    for (int i = 0; i < LANES; ++i) r[i] = __builtin_sqrt(x[i]);
    return r;
}
#endif

static inline bool anyTrue(vi m) {
    long long acc = 0;
    for (int i = 0; i < LANES; ++i) acc |= m[i];
    return acc != 0;
}
static inline bool allTrue(vi m) {
    long long acc = ~0LL;
    for (int i = 0; i < LANES; ++i) acc &= m[i];
    return acc != 0;
}

// Round to nearest without an intrinsic: adding and subtracting 2^52 + 2^51
// forces the fractional bits out under round-to-nearest-even.  Valid for
// |x| < 2^51, which the quadrant index here never approaches.
static inline vd vround(vd x) {
    const vd magic = vsplat(6755399441055744.0);
    vd t = (x + magic) - magic;
    // The trick collapses toward zero for huge inputs; fall back for those.
    return vsel(vabs(x) < vsplat(4.503599627370496e15), t, x);
}

// ---------------------------------------------------------------------------
// sin and cos together
// ---------------------------------------------------------------------------
//
// Cody-Waite reduction to |r| <= pi/4 against a three-part pi/2, then Taylor
// series in r.  Terms are carried far enough that the truncation error is
// ~1e-16 relative at the fold edge, so replacing libm here does not measurably
// change a trajectory -- verified against std::sin/std::cos in --check.
static inline void vsincos(vd x, vd& sinOut, vd& cosOut) {
    const vd TWO_OVER_PI = vsplat(0.636619772367581343076);
    vd kf = vround(x * TWO_OVER_PI);

    // pi/2 split so that kf * each part is exact in double.
    const vd PIO2_1 = vsplat(1.57079632673412561417e+00);
    const vd PIO2_2 = vsplat(6.07710050650619224932e-11);
    const vd PIO2_3 = vsplat(2.02226624879595063154e-21);
    vd r = x - kf * PIO2_1;
    r = r - kf * PIO2_2;
    r = r - kf * PIO2_3;

    vd r2 = r * r;

    // sin(r)/r, Taylor through r^15
    vd s = vsplat(-7.64716373181981647590e-13);         // -1/15!
    s = s * r2 + vsplat( 1.60590438368216145993e-10);   //  1/13!
    s = s * r2 + vsplat(-2.50521083854417187751e-08);   // -1/11!
    s = s * r2 + vsplat( 2.75573192239858906526e-06);   //  1/9!
    s = s * r2 + vsplat(-1.98412698412698412698e-04);   // -1/7!
    s = s * r2 + vsplat( 8.33333333333333333333e-03);   //  1/5!
    s = s * r2 + vsplat(-1.66666666666666666667e-01);   // -1/3!
    vd S = r + r * r2 * s;

    // cos(r), Taylor through r^14
    vd c = vsplat(-1.14707455977297247139e-11);         // -1/14!
    c = c * r2 + vsplat( 2.08767569878680989792e-09);   //  1/12!
    c = c * r2 + vsplat(-2.75573192239858906526e-07);   // -1/10!
    c = c * r2 + vsplat( 2.48015873015873015873e-05);   //  1/8!
    c = c * r2 + vsplat(-1.38888888888888888889e-03);   // -1/6!
    c = c * r2 + vsplat( 4.16666666666666666667e-02);   //  1/4!
    c = c * r2 + vsplat(-5.00000000000000000000e-01);   // -1/2!
    vd C = vsplat(1.0) + r2 * c;

    // Quadrant fold.  k odd swaps the two series; the sign of each follows
    // bit 1 of k (and of k+1 for cosine).
    vi k = __builtin_convertvector(kf, vi);
    vi swap   = (k & isplat(1)) != isplat(0);
    vi sinNeg = (k & isplat(2)) != isplat(0);
    vi cosNeg = ((k + isplat(1)) & isplat(2)) != isplat(0);

    vd sa = vsel(swap, C, S);
    vd ca = vsel(swap, S, C);

    const vi sign = isplat(int64_t(0x8000000000000000ULL));
    sinOut = vbits(ibits(sa) ^ (sinNeg & sign));
    cosOut = vbits(ibits(ca) ^ (cosNeg & sign));
}

// x^(-1/5) for the step controller, lanewise.  Same seed-plus-Newton scheme as
// the scalar path: exponent arithmetic for a rough guess, then three Newton
// iterations on f(y) = y^-5 - x.
static inline vd vInvFifthRoot(vd x) {
    x = vmin(vmax(x, vsplat(1e-30)), vsplat(1e30));
    const vi B = isplat(4607182418800017408LL);          // 1023 << 52
    vi bits = B - (ibits(x) - B) / isplat(5);
    vd y = vbits(bits);
    for (int i = 0; i < 3; ++i) {
        vd y2 = y * y;
        vd y5 = y2 * y2 * y;
        y = y * ((vsplat(6.0) - x * y5) * vsplat(0.2));
    }
    return y;
}
