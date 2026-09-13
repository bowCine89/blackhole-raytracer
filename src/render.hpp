// render.hpp -- the tracing core, shared by the batch renderer (main.cpp) and
// the interactive viewer (viewer.cpp).
//
// Everything here is pure: given a scene and a camera it produces radiance.
// Neither the CLI nor the SDL front end has its own copy of the physics.
#pragma once
#include "core.hpp"
#include "kerr.hpp"
#include "spectrum.hpp"
#include "scene.hpp"

// ---------------------------------------------------------------------------
// Geodesic propagation
// ---------------------------------------------------------------------------
enum class Term { Captured, Escaped, Disk, Stalled };

struct Propagator {
    const Kerr*  kerr;
    const Disk*  disk;
    Real rEscape = 2000;
    Real rCapture = 2;
    Real rtol = 1e-6;      // the render default; --rtol 1e-9 for reference quality
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

    // A transparent outer disk can be crossed several times by a single
    // strongly lensed ray; this only bounds pathological cases.
    constexpr int MAX_CROSSINGS = 32;
    int bounce = 0, crossings = 0;

    for (;;) {
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

        if (++crossings > MAX_CROSSINGS) break;

        // Radiative transfer through the layer.  For an isothermal slab,
        //     I_out = I_in e^-tau + B(T) (1 - e^-tau),
        // so with probability e^-tau the photon passes through untouched and
        // otherwise the layer emits its blackbody.  Sampling it this way is
        // unbiased and collapses to the old opaque disk as tau -> infinity.
        //
        // e_2 is exactly d/dtheta normalised, so the fluid-frame direction
        // cosine against the disk normal is just p_theta / (r nu).  A ray that
        // skims the disk crosses more material: tau_eff = tau_perp / |mu|.
        Real mu = clampf(std::fabs(pth) / (r * nu), 1e-4, 1.0);
        Real trans = std::exp(-disk.opticalDepth(r) / mu);

        if (trans > 0 && rng.uniform() < trans) {
            // Passes through.  Same photon, same momentum, so `shift` is
            // untouched -- only nudge past the plane so the crossing search
            // does not re-fire on the point we just left.
            g.r   = r;
            g.th  = HALF_PI + (pth > 0 ? 1e-11 : -1e-11);
            g.ph  = ph;
            g.pr  = yEnd.y[3];
            g.pth = pth;
            continue;
        }

        shift *= nu;

        Real gs = 1 / shift;
        Real gs2 = gs * gs, gs5 = gs2 * gs2 * gs;
        Real T = disk.temperature(r, ph);
        if (T > 0) radiance += throughput * gs5 * spec::planck(lambda0 * gs, T);

        if (bounce++ >= maxBounces) break;
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
// Where the camera is and where it is pointing.  Kept separate from the
// Camera itself so the interactive viewer can mutate a handful of numbers and
// rebuild the frame, without knowing anything about tetrads.
struct CameraParams {
    Real camR      = 40.0;    // Boyer-Lindquist radius, in M
    Real incDeg    = 80.0;    // inclination from the spin axis; 90 = edge-on
    Real phiDeg    = 0.0;     // azimuth around the spin axis
    Real fovDeg    = 40.0;    // horizontal field of view
    Real yawDeg    = 0.0;     // aim offset, left/right
    Real pitchDeg  = 0.0;     // aim offset, up/down
    Real aspect    = 16.0 / 9.0;
};

struct Camera {
    Real r, th, ph;
    Real tanHalf, aspect;
    Vec3 fwd, up, right;             // in ZAMO frame components (n_r, n_th, n_ph)
    Zamo zamo;

    Camera(const Kerr& k, const CameraParams& p)
        : r(p.camR), th(p.incDeg * PI / 180), ph(p.phiDeg * PI / 180),
          tanHalf(std::tan(0.5 * p.fovDeg * PI / 180)),
          aspect(p.aspect),
          zamo(k, p.camR, p.incDeg * PI / 180)
    {
        fwd   = {-1, 0, 0};          // toward the hole
        up    = {0, -1, 0};          // toward decreasing theta, i.e. "north"
        right = {0, 0, 1};           // direction of increasing phi
        rotate(up,    p.yawDeg   * PI / 180, fwd, right);
        rotate(right, p.pitchDeg * PI / 180, fwd, up);
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

