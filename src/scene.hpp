// scene.hpp -- what the light actually comes from: a relativistic thin
// accretion disk and a background star field.
#pragma once
#include "core.hpp"
#include "kerr.hpp"
#include "spectrum.hpp"

// ---------------------------------------------------------------------------
// Value noise (used only for optional cosmetic structure in the disk)
// ---------------------------------------------------------------------------
inline Real smoothstep5(Real t) { return t * t * t * (t * (t * 6 - 15) + 10); }

// 2-D value noise that is exactly periodic in y with period `periodY` lattice
// cells, so the disk has no seam where phi wraps around.
inline Real valueNoise2(Real x, Real y, int periodY, uint32_t seed) {
    int xi = int(std::floor(x)), yi = int(std::floor(y));
    Real fx = x - xi, fy = y - yi;
    auto at = [&](int i, int j) {
        int jm = ((j % periodY) + periodY) % periodY;
        return hashFloat(hashCombine(uint32_t(i + 8192), uint32_t(jm), seed));
    };
    Real u = smoothstep5(fx), v = smoothstep5(fy);
    Real a = at(xi, yi),     b = at(xi + 1, yi);
    Real c = at(xi, yi + 1), d = at(xi + 1, yi + 1);
    return (a + (b - a) * u) + ((c + (d - c) * u) - (a + (b - a) * u)) * v;
}

inline Real fbm2(Real x, Real y, int periodY, uint32_t seed, int octaves) {
    Real sum = 0, amp = 0.5, norm = 0;
    for (int o = 0; o < octaves; ++o) {
        sum  += amp * valueNoise2(x, y, periodY << o, seed + uint32_t(o) * 977u);
        norm += amp;
        x *= 2; y *= 2; amp *= 0.5;
    }
    return sum / norm;
}

// ---------------------------------------------------------------------------
// Novikov-Thorne / Page-Thorne relativistic thin disk
// ---------------------------------------------------------------------------
//
// Geometrically thin, optically thick, in the equatorial plane, with matter on
// prograde Keplerian circular orbits between the ISCO and rOut.  Each annulus
// radiates as a blackbody at the local effective temperature set by the
// Page & Thorne (1974) dissipation profile, which vanishes at the ISCO
// (zero-torque inner boundary) and falls off as r^-3 far out.
struct Disk {
    Kerr   kerr;
    Real   rIn = 6, rOut = 22;
    Real   tPeak = 12000;        // K, peak effective temperature
    Real   albedo = 0.2;         // grey scattering albedo of the disk surface
    Real   turbulence = 0.15;    // cosmetic temperature mottling, 0 disables
    uint32_t noiseSeed = 12345;

    // --- outer edge ---------------------------------------------------------
    //
    // A real thin disk has no outer edge: the dissipation profile falls off as
    // r^-3 and the disk simply continues out to wherever it is fed.  What ends
    // it observationally is the surface density dropping until the disk is no
    // longer optically thick, at which point it both dims and becomes
    // transparent.  The physically meaningful quantity is therefore the
    // vertical optical depth, not a radius.
    //
    //     tau_perp(r) = min(tauMax, exp((rOut - r) / edgeWidth))
    //
    // so rOut is the photosphere edge (tau_perp = 1) rather than a hard stop,
    // and edgeWidth is the e-folding length of the decline.  Setting
    // edgeWidth = 0 restores the old hard-edged disk exactly.
    Real   edgeWidth = 1.2;      // e-folding width of the optical-depth decline
    Real   tauMax = 30.0;        // interior optical depth (see note below)
    Real   rCut = 22;            // radius beyond which tau is negligible

    // tauMax is capped at 30 rather than the 1e4-1e6 a real disk carries.
    // exp(-30) is 1e-13: already perfectly opaque to within floating-point
    // noise, so the cap is numerically indistinguishable from the real value
    // while keeping the taper width interpretable.

    // Page-Thorne auxiliaries
    Real x0 = 0, x1 = 0, x2 = 0, x3 = 0, aStar = 0;
    Real fluxNorm = 1;           // 1 / max flux, so the profile peaks at 1

    Disk() = default;

    void init(const Kerr& k, Real innerR, Real outerR) {
        kerr = k;
        // The flux formula is singular at exactly a = 0 (a removable 0/0); a
        // hair of spin keeps it finite, with a relative error around 1e-7.
        aStar = k.a;
        if (std::fabs(aStar) < 1e-6) aStar = 1e-6;

        Real rIsco = k.isco(true);
        rIn  = std::max(innerR > 0 ? innerR : rIsco, k.horizon() * 1.001);
        rOut = outerR;

        Real ac = (1.0 / 3.0) * std::acos(clampf(aStar, -1, 1));
        x1 =  2 * std::cos(ac - PI / 3);
        x2 =  2 * std::cos(ac + PI / 3);
        x3 = -2 * std::cos(ac);
        x0 = std::sqrt(rIsco);

        // Normalise so the peak of the profile is 1.
        fluxNorm = 1;
        Real fmax = 0;
        const int N = 4096;
        for (int i = 1; i < N; ++i) {
            Real r = rIn + (rOut - rIn) * i / Real(N);
            fmax = std::max(fmax, rawFlux(r));
        }
        fluxNorm = fmax > 0 ? 1 / fmax : 1;

        // Stop tracking the disk once it transmits all but a fraction of a
        // percent; beyond that it is indistinguishable from empty space.
        rCut = (edgeWidth > 0) ? rOut + edgeWidth * std::log(200.0) : rOut;
    }

    // Geometric extent the integrator has to consider.
    bool contains(Real r) const { return r >= rIn && r <= rCut; }

    // Vertical (face-on) optical depth of the disk at radius r.
    Real opticalDepth(Real r) const {
        if (r < rIn || r > rCut) return 0;
        if (edgeWidth <= 0) return tauMax;
        return std::min(tauMax, std::exp((rOut - r) / edgeWidth));
    }

    // Page & Thorne (1974) eq. (15n): energy flux radiated per unit proper
    // area, up to a constant (Mdot, M) prefactor that the normalisation eats.
    Real rawFlux(Real r) const {
        Real x = std::sqrt(r);
        if (x <= x0) return 0;
        Real denom = x * x * x - 3 * x + 2 * aStar;
        if (denom <= 0) return 0;

        auto term = [&](Real xi, Real xj, Real xk) {
            Real num = 3 * (xi - aStar) * (xi - aStar);
            Real den = xi * (xi - xj) * (xi - xk);
            Real arg = (x - xi) / (x0 - xi);
            if (den == 0 || arg <= 0) return Real(0);
            return num / den * std::log(arg);
        };

        Real bracket = x - x0 - 1.5 * aStar * std::log(x / x0)
                     - term(x1, x2, x3) - term(x2, x1, x3) - term(x3, x1, x2);
        Real f = bracket / (x * x * x * x * denom);
        return std::max<Real>(0, f) * fluxNorm;
    }

    // Local effective temperature.  F = sigma T^4, so T scales as F^(1/4).
    Real temperature(Real r, Real phi) const {
        Real f = rawFlux(r);
        if (f <= 0) return 0;
        Real T = tPeak * std::pow(f, 0.25);
        if (turbulence > 0) {
            // Mottling that shears with radius, as a differentially rotating
            // flow would.  Purely cosmetic -- set --turbulence 0 for the clean
            // Novikov-Thorne profile.
            Real u = std::log(r) * 6.0;
            Real v = (phi + 2.5 * std::pow(r, -1.5) * 40.0) * (16 / TWO_PI);
            Real n = fbm2(u, v, 16, noiseSeed, 4);
            T *= (1 + turbulence * (2 * n - 1));
        }
        return T;
    }

    // Four-velocity of the prograde Keplerian circular orbit at radius r in
    // the equatorial plane, in Boyer-Lindquist components.
    Vec4 orbitVelocity(Real r) const {
        Real a = kerr.a;
        Real r32 = r * std::sqrt(r);
        Real omega = 1 / (r32 + a);
        Real denom = r * r * r - 3 * r * r + 2 * a * r32;
        Real ut = (denom > 0) ? (r32 + a) / std::sqrt(denom) : 1e6;
        Vec4 u;
        u[0] = ut; u[1] = 0; u[2] = 0; u[3] = omega * ut;
        return u;
    }

    Real orbitOmega(Real r) const { return 1 / (r * std::sqrt(r) + kerr.a); }
};

// ---------------------------------------------------------------------------
// Procedural star field
// ---------------------------------------------------------------------------
//
// Stars live in cube-face cells and are generated on demand from an integer
// hash, so the "catalogue" costs no memory and is identical on every thread.
// Each star is a blackbody, so it redshifts through the lens like everything
// else.
struct Sky {
    Real gain       = 1.0;
    Real density    = 0.38;      // fraction of cells holding a star
    Real sigma      = 5.0e-4;    // angular radius (radians)
    int  cells      = 420;       // cells per cube-face edge
    Real bandGain   = 0.10;      // faint galactic-plane glow, 0 disables
    uint32_t seed   = 7777;
    bool enabled    = true;

    // Map a unit direction to (face, u, v) with u,v in [-1,1].
    static void toCube(const Vec3& d, int& face, Real& u, Real& v) {
        Real ax = std::fabs(d.x), ay = std::fabs(d.y), az = std::fabs(d.z);
        if (ax >= ay && ax >= az)      { face = d.x > 0 ? 0 : 1; u = d.y / ax; v = d.z / ax; }
        else if (ay >= az)             { face = d.y > 0 ? 2 : 3; u = d.z / ay; v = d.x / ay; }
        else                           { face = d.z > 0 ? 4 : 5; u = d.x / az; v = d.y / az; }
    }
    // Exact inverse of toCube.  The negative faces must NOT negate u and v:
    // toCube divides by |component|, so (-1, y, z) yields (u, v) = (y, z) and
    // the direction to rebuild is (-1, u, v).  Negating here puts every star on
    // a negative face in the antipodal direction, where the angular cutoff
    // discards it -- silently emptying half the sky.
    static Vec3 fromCube(int face, Real u, Real v) {
        switch (face) {
            case 0: return normalize({ 1,  u,  v});
            case 1: return normalize({-1,  u,  v});
            case 2: return normalize({ v,  1,  u});
            case 3: return normalize({ v, -1,  u});
            case 4: return normalize({ u,  v,  1});
            default:return normalize({ u,  v, -1});
        }
    }

    // Spectral radiance of the background along direction d at wavelength l.
    //
    // `withStars` is false for rays that have already scattered off the disk.
    // Stars are very nearly delta functions, so sampling them through a diffuse
    // bounce produces severe fireflies -- and the signal being sampled is
    // physically negligible anyway: an annulus radiating as a 9000 K blackbody
    // outshines the starlight it reflects by some ten orders of magnitude.
    // Dropping them there removes the variance without removing anything real.
    Real radiance(const Vec3& d, Real lambda, bool withStars = true) const {
        if (!enabled) return 0;
        Real total = 0;

        if (withStars) {
            int face; Real u, v;
            toCube(d, face, u, v);
            Real fu = (u * 0.5 + 0.5) * cells, fv = (v * 0.5 + 0.5) * cells;
            int ci = int(std::floor(fu)), cj = int(std::floor(fv));

            const Real inv2s2 = 1 / (2 * sigma * sigma);
            for (int dj = -1; dj <= 1; ++dj) {
                for (int di = -1; di <= 1; ++di) {
                    int i = ci + di, j = cj + dj;
                    if (i < 0 || j < 0 || i >= cells || j >= cells) continue;
                    uint32_t h = hashCombine(uint32_t(i), uint32_t(j), seed + uint32_t(face) * 9176u);
                    if (hashFloat(h) > density) continue;

                    uint32_t h2 = hashU32(h ^ 0x51ed270bU);
                    uint32_t h3 = hashU32(h2 ^ 0x2545f491U);
                    uint32_t h4 = hashU32(h3 ^ 0x9e3779b9U);

                    Real su = ((i + hashFloat(h2)) / cells) * 2 - 1;
                    Real sv = ((j + hashFloat(h3)) / cells) * 2 - 1;
                    Vec3 sd = fromCube(face, su, sv);

                    // Chord distance approximates the angular separation to well
                    // under a part in 10^6 at these scales.
                    Vec3 diff = d - sd;
                    Real d2 = dot(diff, diff);
                    Real w = d2 * inv2s2;
                    if (w > 12) continue;

                    // Steep flux distribution: many faint stars, few bright ones.
                    Real q = hashFloat(h4);
                    Real flux = q * q * q * q;
                    // Colour temperature, biased toward cool stars with a blue tail.
                    Real tq = hashFloat(hashU32(h4 ^ 0x68bc21ebU));
                    Real T = 2800 + 9000 * tq * tq + (tq > 0.97 ? 18000 * (tq - 0.97) / 0.03 : 0);

                    total += flux * std::exp(-w) * spec::planck(lambda, T);
                }
            }
            total *= gain;
        }

        if (bandGain > 0) {
            // A soft luminous band to stand in for the galactic plane; it
            // gives the lensed sky some large-scale structure to distort.
            Real b = d.y;                                  // band normal is +y
            Real prof = std::exp(-(b * b) / (2 * 0.13 * 0.13));
            Real n = fbm2(d.x * 3.0 + 10, d.z * 3.0 + 10, 64, seed ^ 0xabcdu, 4);
            total += bandGain * prof * (0.35 + 0.65 * n) * spec::planck(lambda, 4200) * 0.02;
        }
        return total;
    }
};
