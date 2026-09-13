// spectrum.hpp -- physically-based colour: Planck emission, CIE 1931 matching
// functions, XYZ -> sRGB, and tone mapping.
//
// The renderer is spectral rather than RGB.  Each camera sample carries a
// single wavelength; the relativistic frequency shift then moves energy
// between wavelengths exactly as it does in reality, so the approaching side
// of the disk really does go blue and the receding side really does go red,
// with the correct hue rather than a hand-tuned tint.
#pragma once
#include "core.hpp"

namespace spec {

constexpr Real LAMBDA_MIN = 380e-9;   // m
constexpr Real LAMBDA_MAX = 780e-9;
constexpr Real LAMBDA_SPAN = LAMBDA_MAX - LAMBDA_MIN;

constexpr Real H_PLANCK = 6.62607015e-34;   // J s
constexpr Real C_LIGHT  = 2.99792458e8;     // m/s
constexpr Real K_BOLTZ  = 1.380649e-23;     // J/K
constexpr Real SIGMA_SB = 5.670374419e-8;   // W m^-2 K^-4

// Spectral radiance of a blackbody per unit wavelength, W m^-2 sr^-1 m^-1.
// Scaled by 1e-13 to keep the numbers in a comfortable floating-point range;
// this is a global constant factor and is absorbed by the exposure control.
inline Real planck(Real lambda, Real T) {
    if (T <= 0) return 0;
    const Real c1 = 2 * H_PLANCK * C_LIGHT * C_LIGHT;
    const Real c2 = H_PLANCK * C_LIGHT / K_BOLTZ;
    Real l5 = lambda * lambda * lambda * lambda * lambda;
    Real x  = c2 / (lambda * T);
    if (x > 700) return 0;                       // exp would overflow; radiance is nil
    return (c1 / l5) / (std::exp(x) - 1) * 1e-13;
}

// CIE 1931 colour matching functions, multi-lobe Gaussian fit of
// Wyman, Sloan & Shirley (JCGT 2013).  Max error ~1%, and it costs a handful
// of exp() calls instead of a table lookup and interpolation.
inline Real gaussLobe(Real x, Real mu, Real s1, Real s2) {
    Real t = (x - mu) * (x < mu ? 1 / s1 : 1 / s2);
    return std::exp(-0.5 * t * t);
}
inline Vec3 cieXYZ(Real lambdaMetres) {
    Real l = lambdaMetres * 1e9;   // the fit is in nanometres
    Real X = 1.056 * gaussLobe(l, 599.8, 37.9, 31.0)
           + 0.362 * gaussLobe(l, 442.0, 16.0, 26.7)
           - 0.065 * gaussLobe(l, 501.1, 20.4, 26.2);
    Real Y = 0.821 * gaussLobe(l, 568.8, 46.9, 40.5)
           + 0.286 * gaussLobe(l, 530.9, 16.3, 31.1);
    Real Z = 1.217 * gaussLobe(l, 437.0, 11.8, 36.0)
           + 0.681 * gaussLobe(l, 459.0, 26.0, 13.8);
    return {X, Y, Z};
}

// Integral of ybar over the sampled band, used to normalise exposure so that
// results are independent of the chosen wavelength range.
inline Real ybarIntegral() {
    Real s = 0;
    const int N = 4001;
    for (int i = 0; i < N; ++i) {
        Real l = LAMBDA_MIN + LAMBDA_SPAN * (i + 0.5) / N;
        s += cieXYZ(l).y;
    }
    return s * (LAMBDA_SPAN / N);
}

// CIE XYZ (D65) -> linear sRGB.
inline Vec3 xyzToLinearSrgb(const Vec3& c) {
    return { 3.2404542 * c.x - 1.5371385 * c.y - 0.4985314 * c.z,
            -0.9692660 * c.x + 1.8760108 * c.y + 0.0415560 * c.z,
             0.0556434 * c.x - 0.2040259 * c.y + 1.0572252 * c.z };
}

// Narkowicz's fitted ACES filmic curve; keeps the enormous dynamic range of
// the inner disk from clipping to a flat white disc.
inline Real acesFilmic(Real x) {
    x = std::max<Real>(0, x);
    return clampf((x * (2.51 * x + 0.03)) / (x * (2.43 * x + 0.59) + 0.14), 0, 1);
}

inline Real srgbEncode(Real u) {
    u = clampf(u, 0, 1);
    return u <= 0.0031308 ? 12.92 * u : 1.055 * std::pow(u, 1 / 2.4) - 0.055;
}

// Pull a strongly saturated colour toward white by the given amount.  Real
// cameras and eyes desaturate at high luminance; without this the hottest
// parts of the disk stay a lurid cyan instead of going white-hot.
inline Vec3 desaturateHighlights(const Vec3& rgb, Real strength) {
    Real m = std::max({rgb.x, rgb.y, rgb.z});
    if (m <= 1 || strength <= 0) return rgb;
    Real t = clampf(strength * (1 - 1 / m), 0, 1);
    return { rgb.x + (m - rgb.x) * t, rgb.y + (m - rgb.y) * t, rgb.z + (m - rgb.z) * t };
}

} // namespace spec
