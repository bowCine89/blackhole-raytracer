// main.cpp -- a spectral Monte Carlo path tracer for the Kerr spacetime.
//
// Rays are null geodesics of a spinning (Kerr) black hole, integrated backward
// from a ZAMO camera.  They can terminate on the horizon, escape to the star
// field, or strike a Novikov-Thorne accretion disk, where they are re-emitted
// as a blackbody and optionally scattered onward -- which is what produces
// "returning radiation", light that leaves the disk, loops around the hole and
// lands back on it.
//
// Radiative transfer is done on the Lorentz invariant I_lambda * lambda^5,
// so gravitational redshift, frame dragging and orbital Doppler beaming all
// fall out of one bookkeeping variable instead of being applied as effects.
#include "core.hpp"
#include "kerr.hpp"
#include "spectrum.hpp"
#include "scene.hpp"
#include "image.hpp"

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <vector>
#include <algorithm>

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------
struct Config {
    int   width = 1280, height = 720;
    int   spp = 128;
    int   maxBounces = 2;
    int   threads = 0;               // 0 = hardware concurrency
    uint64_t seed = 20260913;

    Real  spin = 0.94;               // a/M
    Real  camR = 40.0;
    Real  camIncDeg = 80.0;          // 90 = edge-on
    Real  camPhiDeg = 0.0;
    Real  fovDeg = 40.0;
    Real  yawDeg = 0.0, pitchDeg = 0.0;

    Real  diskIn = -1;               // <0 -> ISCO
    Real  diskOut = 18.0;
    Real  tPeak = 9000.0;            // K
    Real  albedo = 0.20;
    Real  turbulence = 0.15;

    Real  skyGain = 0.20;
    Real  starDensity = 0.16;
    Real  bandGain = 0.15;
    bool  stars = true;

    Real  exposure = 1.0;
    Real  key = 1.3;
    Real  percentile = 0.995;
    Real  bloom = 0.16;
    Real  bloomThreshold = 1.0;
    Real  desat = 0.85;

    Real  rtol = 1e-6;
    int   maxSteps = 60000;

    std::string out = "blackhole.png";
    bool  writePfm = false;
    bool  check = false;
    bool  stats = false;
    bool  quiet = false;
};

// ---------------------------------------------------------------------------
// Geodesic propagation
// ---------------------------------------------------------------------------
enum class Term { Captured, Escaped, Disk, Stalled };

struct Propagator {
    const Kerr*  kerr;
    const Disk*  disk;
    Real rEscape = 2000;
    Real rCapture = 2;
    Real rtol = 1e-9;
    int  maxSteps = 60000;
    bool diskActive = true;

    // Advance the geodesic until it terminates.  On a disk hit, yEnd holds the
    // state interpolated exactly onto the equatorial plane.
    Term run(const Geodesic& g, State& yEnd, uint64_t* stepsOut = nullptr) const {
        DormandPrince dp;
        dp.rtol = rtol;

        State y{{g.r, g.th, g.ph, g.pr, g.pth}};
        RayConst rc(kerr->a, g.E, g.L);
        geodesicRhs(rc, y, dp.k1);

        // Step caps keep each step geometrically sane: no more than ~30% of
        // the current radius, ~0.25 rad in theta, ~0.5 rad in phi.  Without
        // them a very accurate step could still leap across the whole disk.
        auto capFor = [&](const State& s, const State& k1) {
            Real c = 8.0;
            c = std::min(c, 0.30 * std::max<Real>(s.y[0], 1.0) / std::max(std::fabs(k1.y[0]), 1e-30));
            c = std::min(c, 0.25 / std::max(std::fabs(k1.y[1]), 1e-30));
            c = std::min(c, 0.50 / std::max(std::fabs(k1.y[2]), 1e-30));
            return c;
        };

        Real h = capFor(y, dp.k1);
        int rejects = 0;

        int step = 0;
        struct StepTally { uint64_t* out; int* n; ~StepTally() { if (out) *out += uint64_t(*n); } } tally{stepsOut, &step};
        for (; step < maxSteps; ++step) {
            Real cap = capFor(y, dp.k1);
            if (h > cap) h = cap;
            if (h < dp.hMin) h = dp.hMin;

            State yn;
            Real err = dp.trial(rc, y, h, yn);

            if (!(err <= 1.0) && !std::isfinite(err)) {
                h *= 0.25;
                if (++rejects > 40) return Term::Captured;   // hopeless; treat as lost
                continue;
            }
            if (err > 1.0) {
                Real hNew = DormandPrince::nextStep(h, err, dp.hMin, cap);
                if (++rejects > 40 || hNew <= dp.hMin * 1.001) { /* force progress */ }
                else { h = hNew; continue; }
            }
            rejects = 0;

            // Equatorial-plane crossing.  The theta motion has turning points
            // strictly inside (0, pi), so theta never leaves that interval and
            // sign(cos theta) is simply sign(pi/2 - theta).  Testing theta
            // directly keeps the whole crossing search free of trigonometry,
            // which is worth real time: the bisection below would otherwise
            // cost dozens of cosines per crossing.
            if (diskActive) {
                Real c0 = HALF_PI - y.y[1], c1 = HALF_PI - yn.y[1];
                if (c0 * c1 < 0) {
                    dp.buildInterpolant(y, yn, h);
                    Real lo = 0, hi = 1;
                    State ym;
                    // 40 halvings resolve the crossing to ~1e-12 of a step,
                    // far below the accuracy of the step itself.
                    for (int it = 0; it < 40; ++it) {
                        Real mid = 0.5 * (lo + hi);
                        dp.interpolate(mid, ym);
                        if ((HALF_PI - ym.y[1]) * c0 > 0) lo = mid; else hi = mid;
                    }
                    dp.interpolate(0.5 * (lo + hi), ym);
                    if (disk->contains(ym.y[0])) { yEnd = ym; return Term::Disk; }
                }
            }

            y = yn;
            dp.k1 = dp.k7;                                   // FSAL
            h = DormandPrince::nextStep(h, err, dp.hMin, cap);

            if (y.y[0] > rEscape)                  { yEnd = y; return Term::Escaped; }
            if (y.y[3] < 0 && y.y[0] < rCapture)   { yEnd = y; return Term::Captured; }
            if (!std::isfinite(y.y[0]))            { yEnd = y; return Term::Captured; }
        }
        yEnd = y;
        return Term::Stalled;
    }
};

// Propagation direction at the escape radius, as a Cartesian unit vector with
// +z along the spin axis.  Residual bending beyond r = 2000 M is ~1e-5 rad,
// well under one star radius, so this is the asymptotic direction for our
// purposes.
inline Vec3 asymptoticDirection(const Kerr& k, const State& y, Real E, Real L) {
    Real r = y.y[0], th = y.y[1], ph = y.y[2];
    Zamo z(k, r, th);
    Real nu  = (E - z.omega * L) / z.alpha;
    Real nr  = y.y[3] / z.sqrtSoverD / nu;
    Real nth = y.y[4] / z.sqrtS / nu;
    Real nph = L / z.sSqrtAoverS / nu;
    Real s = std::sin(th), c = std::cos(th), cp = std::cos(ph), sp = std::sin(ph);
    Vec3 er{s * cp, s * sp, c}, et{c * cp, c * sp, -s}, ep{-sp, cp, 0};
    return normalize(er * nr + et * nth + ep * nph);
}

// ---------------------------------------------------------------------------
// The path tracer
// ---------------------------------------------------------------------------
//
// Returns the spectral radiance the camera receives at wavelength lambda0.
//
// `shift` tracks nu_local / nu_camera along the path.  Its reciprocal is the
// familiar g-factor, and because I_lambda * lambda^5 is invariant, a source of
// local temperature T contributes g^5 * B_lambda(g * lambda0, T).
inline Real tracePath(const Kerr& kerr, const Disk& disk, const Sky& sky,
                      const Propagator& prop, Geodesic g, Real lambda0,
                      Rng& rng, int maxBounces, uint64_t* steps = nullptr)
{
    Real radiance = 0;
    Real shift = 1;
    // Stays at unity while the albedo is grey and scattering is chosen by
    // Russian roulette (survivors carry full weight).  It is carried explicitly
    // so a wavelength-dependent or non-Lambertian disk BRDF only has to touch
    // the scattering block.
    Real throughput = 1;

    for (int bounce = 0; ; ++bounce) {
        State yEnd;
        Term t = prop.run(g, yEnd, steps);

        if (t == Term::Captured || t == Term::Stalled) break;

        if (t == Term::Escaped) {
            // At infinity the emitter is static, so nu_there = -p.u = E.
            Real nuRatio = shift * g.E;
            if (!(nuRatio > 0)) break;
            Real gs = 1 / nuRatio;
            Real gs2 = gs * gs, gs5 = gs2 * gs2 * gs;
            Vec3 dir = asymptoticDirection(kerr, yEnd, g.E, g.L);
            radiance += throughput * gs5 * sky.radiance(dir, lambda0 * gs, bounce == 0);
            break;
        }

        // --- disk surface -------------------------------------------------
        Real r  = yEnd.y[0];
        Real ph = yEnd.y[2];
        Real pth = yEnd.y[4];

        Vec4 u   = disk.orbitVelocity(r);
        Real om  = disk.orbitOmega(r);
        Real nu  = u[0] * (g.E - om * g.L);          // -p.u in the fluid frame
        if (!(nu > 0)) break;
        shift *= nu;

        Real gs = 1 / shift;
        Real gs2 = gs * gs, gs5 = gs2 * gs2 * gs;
        Real T = disk.temperature(r, ph);
        if (T > 0) radiance += throughput * gs5 * spec::planck(lambda0 * gs, T);

        if (bounce >= maxBounces) break;
        // Russian roulette on the grey albedo: survivors keep unit weight.
        if (rng.uniform() >= disk.albedo) break;

        // Re-emit into the hemisphere the light came from.  pth > 0 means the
        // backward ray was heading south, so the light arrived from the north.
        Real side = (pth > 0) ? -1.0 : 1.0;

        Metric m = kerr.metricAt(r, PI / 2);
        Vec4 e[4];
        buildTetrad(m, u, e);                        // e1 radial, e2 polar, e3 azimuthal

        Vec3 s = cosineHemisphere(rng.uniform(), rng.uniform());
        Vec4 p;                                      // contravariant, unit local energy
        for (int i = 0; i < 4; ++i)
            p[i] = e[0][i] + s.x * e[1][i] + s.y * e[3][i] + (side * s.z) * e[2][i];

        Vec4 pl = m.lower(p);
        g.r   = r;
        g.th  = PI / 2 + side * 1e-11;               // nudge off the plane
        g.ph  = ph;
        g.pr  = pl[1];
        g.pth = pl[2];
        g.E   = -pl[0];
        g.L   = pl[3];
        if (!(g.E > 0)) break;
    }
    return radiance;
}

// ---------------------------------------------------------------------------
// Camera
// ---------------------------------------------------------------------------
struct Camera {
    Real r, th, ph;
    Real tanHalf, aspect;
    Vec3 fwd, up, right;             // in ZAMO frame components (n_r, n_th, n_ph)
    Zamo zamo;

    Camera(const Kerr& k, const Config& c)
        : r(c.camR), th(c.camIncDeg * PI / 180), ph(c.camPhiDeg * PI / 180),
          tanHalf(std::tan(0.5 * c.fovDeg * PI / 180)),
          aspect(Real(c.width) / c.height),
          zamo(k, c.camR, c.camIncDeg * PI / 180)
    {
        fwd   = {-1, 0, 0};          // toward the hole
        up    = {0, -1, 0};          // toward decreasing theta, i.e. "north"
        right = {0, 0, 1};           // direction of increasing phi
        rotate(up,    c.yawDeg   * PI / 180, fwd, right);
        rotate(right, c.pitchDeg * PI / 180, fwd, up);
    }

    // Rotate a and b about the given axis (Rodrigues).
    static void rotate(const Vec3& axis, Real ang, Vec3& a, Vec3& b) {
        if (ang == 0) return;
        Real ca = std::cos(ang), sa = std::sin(ang);
        auto rot = [&](const Vec3& v) {
            return v * ca + cross(axis, v) * sa + axis * (dot(axis, v) * (1 - ca));
        };
        a = normalize(rot(a));
        b = normalize(rot(b));
    }

    Geodesic ray(Real sx, Real sy) const {
        Vec3 d = normalize(fwd + right * (sx * tanHalf * aspect) + up * (sy * tanHalf));
        return zamo.emit(r, th, ph, d);
    }
};

// ---------------------------------------------------------------------------
// Post-processing
// ---------------------------------------------------------------------------
static inline int clampi(int v, int lo, int hi) { return v < lo ? lo : (v > hi ? hi : v); }

static void boxBlurH(std::vector<float>& buf, int w, int h, int rad) {
    std::vector<float> tmp(size_t(w) * h * 3);
    Real inv = 1.0 / (2 * rad + 1);
    for (int y = 0; y < h; ++y) {
        for (int ch = 0; ch < 3; ++ch) {
            double sum = 0;
            auto at = [&](int x) { return double(buf[(size_t(y) * w + clampi(x, 0, w - 1)) * 3 + ch]); };
            for (int x = -rad; x <= rad; ++x) sum += at(x);
            for (int x = 0; x < w; ++x) {
                tmp[(size_t(y) * w + x) * 3 + ch] = float(sum * inv);
                sum += at(x + rad + 1) - at(x - rad);
            }
        }
    }
    buf.swap(tmp);
}
static void boxBlurV(std::vector<float>& buf, int w, int h, int rad) {
    std::vector<float> tmp(size_t(w) * h * 3);
    Real inv = 1.0 / (2 * rad + 1);
    for (int x = 0; x < w; ++x) {
        for (int ch = 0; ch < 3; ++ch) {
            double sum = 0;
            auto at = [&](int y) { return double(buf[(size_t(clampi(y, 0, h - 1)) * w + x) * 3 + ch]); };
            for (int y = -rad; y <= rad; ++y) sum += at(y);
            for (int y = 0; y < h; ++y) {
                tmp[(size_t(y) * w + x) * 3 + ch] = float(sum * inv);
                sum += at(y + rad + 1) - at(y - rad);
            }
        }
    }
    buf.swap(tmp);
}

// ---------------------------------------------------------------------------
// Self-test: validate the geodesic engine against closed-form Kerr results.
// ---------------------------------------------------------------------------
static int runChecks() {
    std::printf("Kerr geodesic engine self-test\n");
    std::printf("------------------------------\n");
    int failures = 0;
    auto report = [&](const char* name, Real got, Real want, Real tol) {
        Real err = std::fabs(got - want) / std::max<Real>(std::fabs(want), 1e-12);
        bool ok = err <= tol;
        if (!ok) ++failures;
        std::printf("  %-46s %14.9f  expected %14.9f  relerr %.2e  %s\n",
                    name, got, want, err, ok ? "ok" : "FAIL");
    };

    // --- horizons, photon orbits, ISCO --------------------------------------
    report("ISCO   a=0     (prograde)",   Kerr(0.0).isco(true),  6.0, 1e-12);
    report("ISCO   a=0.9   (prograde)",   Kerr(0.9).isco(true),  2.320883043, 1e-8);
    report("ISCO   a=0.9   (retrograde)", Kerr(0.9).isco(false), 8.717352279, 1e-8);
    report("photon a=0     (prograde)",   Kerr(0.0).photonOrbit(true), 3.0, 1e-12);
    // Spin is clamped just below extremal, so r_ph(a=1) is not reachable
    // exactly.  Check instead that the closed form satisfies its own defining
    // equation r^(3/2) - 3 r^(1/2) +- 2a = 0, which is an independent test.
    {
        Real worst = 0;
        for (Real a : {0.0, 0.3, 0.7, 0.9, 0.99, 0.9999, -0.5, -0.95}) {
            Kerr k(a);
            for (bool pro : {true, false}) {
                Real r = k.photonOrbit(pro);
                Real res = r * std::sqrt(r) - 3 * std::sqrt(r) + (pro ? 2 * k.a : -2 * k.a);
                worst = std::max(worst, std::fabs(res));
            }
        }
        bool ok = worst < 1e-12;
        if (!ok) ++failures;
        std::printf("  %-46s %14.3e  (tolerance %9.1e)              %s\n",
                    "photon orbit defining-equation residual", worst, 1e-12, ok ? "ok" : "FAIL");
    }

    // --- critical impact parameter by shooting equatorial rays --------------
    // A photon aimed from r = 1000 with impact parameter b either falls in or
    // turns around.  Bisecting on that boundary recovers b_crit, which has a
    // closed form:  b = -(r^3 - 3r^2 + a^2 r + a^2) / (a (r - 1))  at r = r_ph.
    auto criticalB = [](Real a, bool prograde) {
        Kerr k(a);
        Disk none;                       // empty: rIn > rOut so nothing is hit
        none.rIn = 1e30; none.rOut = -1e30;
        Propagator prop;
        prop.kerr = &k; prop.disk = &none;
        prop.rEscape = 2000;
        prop.rCapture = std::max(k.horizon() * 1.0005, k.photonOrbit(true));
        prop.rtol = 1e-11;
        prop.maxSteps = 400000;
        prop.diskActive = false;

        auto captured = [&](Real b) {
            Real r0 = 1000;
            Real D = k.delta(r0);
            Real P = (r0 * r0 + a * a) - a * b;
            Real pr2 = (P * P / D - (b - a) * (b - a)) / D;
            if (pr2 <= 0) return false;
            Geodesic g;
            g.r = r0; g.th = PI / 2; g.ph = 0; g.E = 1; g.L = b;
            g.pth = 0; g.pr = -std::sqrt(pr2);
            State y;
            return prop.run(g, y) != Term::Escaped;
        };
        // Captured for |b| below critical, escapes above.
        Real lo = prograde ? 0.5 : -0.5, hi = prograde ? 12.0 : -12.0;
        for (int i = 0; i < 60; ++i) {
            Real mid = 0.5 * (lo + hi);
            if (captured(mid)) lo = mid; else hi = mid;
        }
        return 0.5 * (lo + hi);
    };
    auto analyticB = [](Real a, bool prograde) {
        Kerr k(a);
        Real r = k.photonOrbit(prograde);
        if (std::fabs(a) < 1e-9) return 3 * std::sqrt(3.0) * (prograde ? 1 : -1);
        return -(r * r * r - 3 * r * r + a * a * r + a * a) / (a * (r - 1));
    };
    report("b_crit a=0                  (3 sqrt 3)", criticalB(1e-9, true), 3 * std::sqrt(3.0), 2e-7);
    report("b_crit a=0.9   (prograde)",   criticalB(0.9, true),   analyticB(0.9, true),   2e-7);
    report("b_crit a=0.9   (retrograde)", criticalB(0.9, false),  analyticB(0.9, false),  2e-7);
    report("b_crit a=0.998 (prograde)",   criticalB(0.998, true), analyticB(0.998, true), 2e-6);

    // --- weak-field light deflection ----------------------------------------
    // End-to-end test of the integrator *and* the asymptotic-direction
    // extraction: a Schwarzschild ray with impact parameter b is deflected by
    //     4M/b + (15 pi/4)(M/b)^2 + (128/3)(M/b)^3 + O((M/b)^4).
    {
        Kerr k(1e-9);
        Disk none; none.rIn = 1e30; none.rOut = -1e30;
        Propagator prop;
        prop.kerr = &k; prop.disk = &none;
        prop.rEscape = 1e6;
        prop.rCapture = std::max(k.horizon() * 1.0005, k.photonOrbit(true));
        prop.rtol = 1e-12; prop.maxSteps = 400000; prop.diskActive = false;

        Real worst = 0;
        for (Real b : {80.0, 200.0, 600.0}) {
            Real r0 = 1e6;
            Real D = k.delta(r0);
            Real P = (r0 * r0 + k.a * k.a) - k.a * b;
            Real pr2 = (P * P / D - (b - k.a) * (b - k.a)) / D;
            Geodesic g;
            g.r = r0; g.th = PI / 2; g.ph = 0; g.E = 1; g.L = b;
            g.pth = 0; g.pr = -std::sqrt(pr2);

            State y0{{g.r, g.th, g.ph, g.pr, g.pth}};
            Vec3 dIn = asymptoticDirection(k, y0, g.E, g.L);
            State y1;
            if (prop.run(g, y1) != Term::Escaped) { ++failures; continue; }
            Vec3 dOut = asymptoticDirection(k, y1, g.E, g.L);

            Real got = std::acos(clampf(dot(dIn, dOut), -1, 1));
            Real m = 1.0 / b;
            Real want = 4 * m + (15 * PI / 4) * m * m + (128.0 / 3.0) * m * m * m;
            worst = std::max(worst, std::fabs(got - want) / want);
        }
        bool ok = worst < 1e-4;
        if (!ok) ++failures;
        std::printf("  %-46s %14.3e  (tolerance %9.1e)              %s\n",
                    "deflection vs post-Newtonian series (rel)", worst, 1e-4, ok ? "ok" : "FAIL");
    }

    // --- the null constraint must stay zero along a traced ray --------------
    {
        Kerr k(0.9);
        Disk none; none.rIn = 1e30; none.rOut = -1e30;
        Propagator prop;
        prop.kerr = &k; prop.disk = &none;
        prop.rEscape = 2000;
        prop.rCapture = std::max(k.horizon() * 1.0005, k.photonOrbit(true));
        prop.rtol = 1e-10; prop.maxSteps = 200000; prop.diskActive = false;

        Config c; Camera cam(k, c);
        Real worst = 0;
        for (int i = 0; i < 400; ++i) {
            Real sx = -1 + 2.0 * ((i * 37) % 400) / 400.0;
            Real sy = -1 + 2.0 * ((i * 91) % 400) / 400.0;
            Geodesic g = cam.ray(sx * 0.9, sy * 0.9);
            State y;
            prop.run(g, y);
            // F is a difference of terms each of order E^2 r^2, so normalise by
            // that: what is being measured is the relative cancellation, and
            // double precision alone costs ~1e-10 of it per evaluation at the
            // escape radius.
            Real scale = std::max<Real>(1, g.E * g.E * (y.y[0] * y.y[0] + k.a * k.a));
            worst = std::max(worst, std::fabs(nullConstraint(k, y, g.E, g.L)) / scale);
        }
        bool ok = worst < 1e-6;
        if (!ok) ++failures;
        std::printf("  %-46s %14.3e  (tolerance %9.1e)              %s\n",
                    "null constraint drift, 400 rays (relative)", worst, 1e-6, ok ? "ok" : "FAIL");
    }

    // --- disk kinematics ----------------------------------------------------
    {
        Kerr k(0.9);
        Disk d; d.init(k, -1, 30);
        // Far from the hole the orbital speed tends to the Newtonian 1/sqrt(r)
        // and u^t -> 1 / sqrt(1 - v^2 - 2/r); check the low-order expansion.
        Real r = 1000;
        Vec4 u = d.orbitVelocity(r);
        Real expected = 1.0 / std::sqrt(1 - 3.0 / r);      // exact for a=0, ok to O(a/r^1.5)
        report("disk u^t at r=1000 (weak field)", u[0], expected, 1e-4);
        Real rIsco = k.isco(true);
        report("disk inner edge == ISCO", d.rIn, rIsco, 1e-12);
        bool zeroTorque = d.rawFlux(rIsco * 1.0000001) < 1e-6;
        if (!zeroTorque) ++failures;
        std::printf("  %-46s %14s                                        %s\n",
                    "Page-Thorne flux vanishes at the ISCO", "", zeroTorque ? "ok" : "FAIL");
    }

    std::printf("------------------------------\n");
    std::printf("%s (%d failure%s)\n", failures ? "FAILED" : "ALL CHECKS PASSED",
                failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------
static void usage() {
    std::printf(
"kerr -- spectral path tracer for a spinning black hole with an accretion disk\n"
"\n"
"Image\n"
"  --width N --height N        output resolution           (1280 720)\n"
"  --spp N                     samples per pixel           (128)\n"
"  --bounces N                 disk scattering bounces     (2)\n"
"  --threads N                 worker threads              (all cores)\n"
"  --seed N                    RNG seed\n"
"  --out FILE                  output .png                 (blackhole.png)\n"
"  --pfm                       also write linear HDR .pfm\n"
"\n"
"Black hole and camera\n"
"  --spin A                    dimensionless spin a/M, |A|<1   (0.94)\n"
"  --dist R                    camera radius in M              (40)\n"
"  --inc DEG                   inclination, 90 = edge-on       (80)\n"
"  --cam-phi DEG               camera azimuth                  (0)\n"
"  --fov DEG                   horizontal field of view        (40)\n"
"  --yaw DEG --pitch DEG       aim offset from the hole        (0 0)\n"
"\n"
"Accretion disk\n"
"  --rin R                     inner radius, <0 means ISCO     (ISCO)\n"
"  --rout R                    outer radius                    (18)\n"
"  --tpeak K                   peak effective temperature      (9000)\n"
"  --albedo A                  disk scattering albedo          (0.20)\n"
"  --turbulence F              cosmetic mottling, 0 = pure NT  (0.15)\n"
"\n"
"Sky\n"
"  --sky-gain F                star brightness                 (0.20)\n"
"  --star-density F            fraction of sky cells with a star (0.16)\n"
"  --band F                    galactic-plane glow             (0.15)\n"
"  --nostars                   black background\n"
"\n"
"Tone mapping\n"
"  --exposure F                exposure multiplier             (1.0)\n"
"  --key F                     target level for the bright end (2.2)\n"
"  --bloom F                   bloom strength, 0 disables      (0.16)\n"
"  --desat F                   highlight desaturation          (0.85)\n"
"\n"
"Accuracy\n"
"  --rtol F                    integrator relative tolerance   (1e-6)\n"
"  --max-steps N               step budget per ray             (60000)\n"
"\n"
"Other\n"
"  --preview                   fast low-resolution settings\n"
"  --check                     run the geodesic self-test and exit\n"
"  --quiet, --help\n");
}

int main(int argc, char** argv) {
    Config c;
    auto needF = [&](int& i) -> Real {
        if (i + 1 >= argc) { std::fprintf(stderr, "missing value for %s\n", argv[i]); std::exit(2); }
        return std::atof(argv[++i]);
    };
    auto needI = [&](int& i) -> int {
        if (i + 1 >= argc) { std::fprintf(stderr, "missing value for %s\n", argv[i]); std::exit(2); }
        return std::atoi(argv[++i]);
    };

    for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];
        if      (s == "--width")       c.width = needI(i);
        else if (s == "--height")      c.height = needI(i);
        else if (s == "--spp")         c.spp = needI(i);
        else if (s == "--bounces")     c.maxBounces = needI(i);
        else if (s == "--threads")     c.threads = needI(i);
        else if (s == "--seed")        c.seed = uint64_t(needI(i));
        else if (s == "--out")         c.out = argv[++i];
        else if (s == "--pfm")         c.writePfm = true;
        else if (s == "--spin")        c.spin = needF(i);
        else if (s == "--dist")        c.camR = needF(i);
        else if (s == "--inc")         c.camIncDeg = needF(i);
        else if (s == "--cam-phi")     c.camPhiDeg = needF(i);
        else if (s == "--fov")         c.fovDeg = needF(i);
        else if (s == "--yaw")         c.yawDeg = needF(i);
        else if (s == "--pitch")       c.pitchDeg = needF(i);
        else if (s == "--rin")         c.diskIn = needF(i);
        else if (s == "--rout")        c.diskOut = needF(i);
        else if (s == "--tpeak")       c.tPeak = needF(i);
        else if (s == "--albedo")      c.albedo = needF(i);
        else if (s == "--turbulence")  c.turbulence = needF(i);
        else if (s == "--sky-gain")    c.skyGain = needF(i);
        else if (s == "--star-density")c.starDensity = needF(i);
        else if (s == "--band")        c.bandGain = needF(i);
        else if (s == "--nostars")     c.stars = false;
        else if (s == "--exposure")    c.exposure = needF(i);
        else if (s == "--key")         c.key = needF(i);
        else if (s == "--bloom")       c.bloom = needF(i);
        else if (s == "--desat")       c.desat = needF(i);
        else if (s == "--rtol")        c.rtol = needF(i);
        else if (s == "--max-steps")   c.maxSteps = needI(i);
        else if (s == "--preview")   { c.width = 480; c.height = 270; c.spp = 24; }
        else if (s == "--check")       c.check = true;
        else if (s == "--quiet")       c.quiet = true;
        else if (s == "--stats")       c.stats = true;
        else if (s == "--help" || s == "-h") { usage(); return 0; }
        else { std::fprintf(stderr, "unknown option: %s\n", s.c_str()); usage(); return 2; }
    }

    if (c.check) return runChecks();

    // --- scene ------------------------------------------------------------
    Kerr kerr(c.spin);
    Disk disk;
    disk.tPeak = c.tPeak;
    disk.albedo = clampf(c.albedo, 0, 0.95);
    disk.turbulence = c.turbulence;
    disk.noiseSeed = uint32_t(c.seed * 2654435761u + 17u);
    disk.init(kerr, c.diskIn, c.diskOut);

    Sky sky;
    sky.gain = c.skyGain;
    sky.density = c.starDensity;
    sky.bandGain = c.bandGain;
    sky.enabled = c.stars;
    sky.seed = uint32_t(c.seed * 40503u + 991u);

    Camera cam(kerr, c);

    // Keep stars roughly a pixel across whatever the resolution is, otherwise
    // a sub-pixel point source turns into chromatic speckle under spectral
    // sampling instead of a clean antialiased dot.
    Real pixelAngle = 2 * cam.tanHalf * cam.aspect / c.width;
    sky.sigma = std::max<Real>(sky.sigma, 0.70 * pixelAngle);

    Propagator prop;
    prop.kerr = &kerr;
    prop.disk = &disk;
    prop.rEscape = std::max<Real>(2000.0, c.camR * 8);
    prop.rCapture = std::max(kerr.horizon() * 1.0005,
                             std::min(kerr.photonOrbit(true), disk.rIn));
    prop.rtol = c.rtol;
    prop.maxSteps = c.maxSteps;

    int nThreads = c.threads > 0 ? c.threads : int(std::thread::hardware_concurrency());
    if (nThreads < 1) nThreads = 1;

    if (!c.quiet) {
        std::printf("Kerr path tracer\n");
        std::printf("  spin a/M          %.4f   (horizon %.4f M, ISCO %.4f M, photon ring %.4f M)\n",
                    kerr.a, kerr.horizon(), kerr.isco(true), kerr.photonOrbit(true));
        std::printf("  camera            r = %.1f M, inclination %.1f deg, fov %.1f deg\n",
                    c.camR, c.camIncDeg, c.fovDeg);
        std::printf("  disk              %.3f M -> %.1f M, T_peak %.0f K, albedo %.2f\n",
                    disk.rIn, disk.rOut, disk.tPeak, disk.albedo);
        std::printf("  image             %d x %d, %d spp, %d bounce%s, %d threads\n",
                    c.width, c.height, c.spp, c.maxBounces, c.maxBounces == 1 ? "" : "s", nThreads);
        std::fflush(stdout);
    }

    // --- render -----------------------------------------------------------
    const int W = c.width, H = c.height;
    std::vector<double> xyzBuf(size_t(W) * H * 3, 0.0);

    const int TILE = 16;
    const int tilesX = (W + TILE - 1) / TILE, tilesY = (H + TILE - 1) / TILE;
    const int nTiles = tilesX * tilesY;
    std::atomic<int> nextTile{0}, doneTiles{0};
    std::atomic<uint64_t> raysTraced{0}, stepsTaken{0};

    const Real ybarInt = spec::ybarIntegral();
    const Real lamScale = spec::LAMBDA_SPAN / ybarInt;

    auto t0 = std::chrono::steady_clock::now();

    auto worker = [&](int tid) {
        uint64_t localRays = 0, localSteps = 0;
        for (;;) {
            int ti = nextTile.fetch_add(1, std::memory_order_relaxed);
            if (ti >= nTiles) break;
            int tx = (ti % tilesX) * TILE, ty = (ti / tilesX) * TILE;
            int x1 = std::min(tx + TILE, W), y1 = std::min(ty + TILE, H);

            for (int y = ty; y < y1; ++y) {
                for (int x = tx; x < x1; ++x) {
                    Rng rng(uint64_t(y) * W + x + 1, c.seed);
                    // Cranley-Patterson rotation of a stratified lattice: the
                    // pixel filter and the wavelength are both well spread.
                    Real rx = rng.uniform(), ry = rng.uniform(), rl = rng.uniform();
                    Vec3 xyz{0, 0, 0};

                    for (int s = 0; s < c.spp; ++s) {
                        Real u1 = rx + 0.7548776662 * (s + 1); u1 -= std::floor(u1);
                        Real u2 = ry + 0.5698402910 * (s + 1); u2 -= std::floor(u2);
                        Real ul = rl + (s + 0.5) / c.spp;      ul -= std::floor(ul);

                        Real lambda = spec::LAMBDA_MIN + spec::LAMBDA_SPAN * ul;
                        Real sx = 2.0 * (x + u1) / W - 1.0;
                        Real sy = 1.0 - 2.0 * (y + u2) / H;

                        Geodesic g = cam.ray(sx, sy);
                        Real L = tracePath(kerr, disk, sky, prop, g, lambda, rng, c.maxBounces, &localSteps);
                        if (L > 0) xyz += spec::cieXYZ(lambda) * L;
                        ++localRays;
                    }

                    Real inv = lamScale / c.spp;
                    size_t o = (size_t(y) * W + x) * 3;
                    xyzBuf[o + 0] = xyz.x * inv;
                    xyzBuf[o + 1] = xyz.y * inv;
                    xyzBuf[o + 2] = xyz.z * inv;
                }
            }
            doneTiles.fetch_add(1, std::memory_order_relaxed);
        }
        raysTraced.fetch_add(localRays, std::memory_order_relaxed);
        stepsTaken.fetch_add(localSteps, std::memory_order_relaxed);
    };

    std::vector<std::thread> pool;
    pool.reserve(nThreads);
    for (int t = 0; t < nThreads; ++t) pool.emplace_back(worker, t);

    if (!c.quiet) {
        int last = -1;
        while (true) {
            int d = doneTiles.load(std::memory_order_relaxed);
            int pct = int(100.0 * d / nTiles);
            if (pct != last) {
                std::printf("\r  rendering        %3d%%", pct);
                std::fflush(stdout);
                last = pct;
            }
            if (d >= nTiles) break;
            std::this_thread::sleep_for(std::chrono::milliseconds(120));
        }
    }
    for (auto& t : pool) t.join();

    auto t1 = std::chrono::steady_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();
    if (!c.quiet) {
        double rays = double(raysTraced.load());
        double steps = double(stepsTaken.load());
        std::printf("\r  rendering        done in %.2f s  (%.2f Mrays/s, %.0f steps/ray, %.1f Msteps/s)\n",
                    secs, rays / secs / 1e6, steps / rays, steps / secs / 1e6);
    }

    // --- exposure ---------------------------------------------------------
    std::vector<double> lum;
    lum.reserve(size_t(W) * H);
    for (size_t i = 0; i < size_t(W) * H; ++i) lum.push_back(xyzBuf[i * 3 + 1]);
    std::sort(lum.begin(), lum.end());
    auto pct = [&](double p) { return lum[std::min(lum.size() - 1, size_t(lum.size() * p))]; };
    double keyLum = pct(c.percentile);
    double scale = (keyLum > 0 ? c.key / keyLum : 1.0) * c.exposure;

    if (c.stats) {
        std::printf("  luminance        p50 %.4g  p90 %.4g  p99 %.4g  p99.9 %.4g  max %.4g\n",
                    pct(0.50), pct(0.90), pct(0.99), pct(0.999), lum.back());
        std::printf("  exposure         anchor p%.2f = %.4g -> key %.3g  (scale %.4g)\n",
                    c.percentile * 100, keyLum, c.key, scale);
    }

    // --- XYZ -> linear sRGB ------------------------------------------------
    std::vector<float> rgb(size_t(W) * H * 3);
    for (size_t i = 0; i < size_t(W) * H; ++i) {
        Vec3 v = spec::xyzToLinearSrgb({xyzBuf[i * 3] * scale,
                                        xyzBuf[i * 3 + 1] * scale,
                                        xyzBuf[i * 3 + 2] * scale});
        rgb[i * 3 + 0] = float(std::max<Real>(0, v.x));
        rgb[i * 3 + 1] = float(std::max<Real>(0, v.y));
        rgb[i * 3 + 2] = float(std::max<Real>(0, v.z));
    }

    if (c.writePfm) img::writePfm(c.out + ".pfm", W, H, rgb);

    // --- bloom -------------------------------------------------------------
    if (c.bloom > 0) {
        std::vector<float> bright(size_t(W) * H * 3);
        for (size_t i = 0; i < size_t(W) * H * 3; ++i)
            bright[i] = float(std::max<Real>(0, rgb[i] - c.bloomThreshold));
        int rad = std::max(2, W / 110);
        for (int pass = 0; pass < 3; ++pass) { boxBlurH(bright, W, H, rad); boxBlurV(bright, W, H, rad); }
        for (size_t i = 0; i < size_t(W) * H * 3; ++i)
            rgb[i] = float(rgb[i] + c.bloom * bright[i]);
    }

    // --- tone map and encode ----------------------------------------------
    std::vector<uint8_t> out8(size_t(W) * H * 3);
    for (size_t i = 0; i < size_t(W) * H; ++i) {
        Vec3 v{rgb[i * 3], rgb[i * 3 + 1], rgb[i * 3 + 2]};
        v = spec::desaturateHighlights(v, c.desat);
        Real r = spec::acesFilmic(v.x), g = spec::acesFilmic(v.y), b = spec::acesFilmic(v.z);
        out8[i * 3 + 0] = uint8_t(clampf(spec::srgbEncode(r) * 255.0 + 0.5, 0, 255));
        out8[i * 3 + 1] = uint8_t(clampf(spec::srgbEncode(g) * 255.0 + 0.5, 0, 255));
        out8[i * 3 + 2] = uint8_t(clampf(spec::srgbEncode(b) * 255.0 + 0.5, 0, 255));
    }

    if (!img::writePng(c.out, W, H, out8)) {
        std::fprintf(stderr, "failed to write %s\n", c.out.c_str());
        return 1;
    }
    if (!c.quiet) std::printf("  wrote            %s\n", c.out.c_str());
    return 0;
}
