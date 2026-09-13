// packet.hpp -- eight rays per thread, one per SIMD lane.
//
// The scalar tracer is latency-bound: a single ray's Dormand-Prince stages are
// strictly sequential, so each right-hand side waits on the previous one's
// sincos and divide and the FPU sits mostly idle.  Nothing about *one* ray can
// fix that.  Eight independent rays in eight lanes can, because their dependency
// chains are unrelated.
//
// Everything is the same mathematics as render.hpp, rewritten so that control
// flow becomes data flow: a lane that wants to reject a step, or has hit the
// disk, or has fallen through the horizon, is handled by a mask instead of a
// branch.  The loop runs until every lane has terminated.
//
// Coherence comes for free here.  A packet is eight samples of the *same pixel*,
// differing only by sub-pixel jitter and wavelength -- and wavelength does not
// enter the geodesic at all, only the emission -- so the eight trajectories are
// nearly identical and lanes terminate at nearly the same step.
#pragma once
#include "simd.hpp"
#include "render.hpp"

struct PacketState { vd y[NSTATE]; };      // r, theta, phi, p_r, p_theta, t

// Per-lane ray constants, the vector twin of RayConst.
struct PacketConst {
    vd a, a2, E, L, aL, aE, L2, a2E2;
};

static inline PacketConst packConst(Real spin, const vd& E, const vd& L) {
    PacketConst c;
    c.a    = vsplat(spin);
    c.a2   = vsplat(spin * spin);
    c.E    = E;
    c.L    = L;
    c.aL   = c.a * L;
    c.aE   = c.a * E;
    c.L2   = L * L;
    c.a2E2 = c.a2 * E * E;
    return c;
}

// Lanewise twin of geodesicRhs.  Same expression, same flop count per lane.
static inline void geodesicRhsV(const PacketConst& c, const PacketState& s, PacketState& out) {
    const vd r = s.y[0], th = s.y[1], pr = s.y[3], pth = s.y[4];

    vd sth, cth;
    vsincos(th, sth, cth);

    // Keep the polar axis out of the denominators, exactly as the scalar path.
    const vd tiny = vsplat(1e-8);
    vi small = vabs(sth) < tiny;
    sth = vsel(small, vsel(sth < vsplat(0.0), vneg(tiny), tiny), sth);

    vd invS  = vsplat(1.0) / sth;
    vd invS2 = invS * invS;

    vd r2    = r * r;
    vd D     = r2 - vsplat(2.0) * r + c.a2;
    vd dD    = vsplat(2.0) * r - vsplat(2.0);
    vd invD  = vsplat(1.0) / D;
    vd PinvD = (c.E * (r2 + c.a2) - c.aL) * invD;

    vd s2 = sth * sth;

    out.y[0] = D * pr;
    out.y[1] = pth;
    out.y[2] = c.L * invS2 - c.aE + c.a * PinvD;
    out.y[3] = vsplat(-0.5) * (dD * pr * pr - vsplat(4.0) * r * c.E * PinvD + PinvD * PinvD * dD);
    out.y[4] = cth * (c.L2 * invS2 * invS - c.a2E2 * sth);
    out.y[5] = c.a * (c.L - c.aE * s2) + (r2 + c.a2) * PinvD;
}

// ---------------------------------------------------------------------------
// Dormand-Prince 5(4) across eight lanes, each with its own step size
// ---------------------------------------------------------------------------
struct PacketDP {
    vd rtol = vsplat(1e-6), atol = vsplat(1e-12);
    vd hMin = vsplat(1e-10);

    PacketState k1{}, k3{}, k4{}, k5{}, k6{}, k7{};
    PacketState c1{}, c2{}, c3{}, c4{}, c5{};

    vd trial(const PacketConst& rc, const PacketState& y, vd h, PacketState& yOut) {
        constexpr double
            a21 = 1.0 / 5,
            a31 = 3.0 / 40,       a32 = 9.0 / 40,
            a41 = 44.0 / 45,      a42 = -56.0 / 15,      a43 = 32.0 / 9,
            a51 = 19372.0 / 6561, a52 = -25360.0 / 2187, a53 = 64448.0 / 6561,  a54 = -212.0 / 729,
            a61 = 9017.0 / 3168,  a62 = -355.0 / 33,     a63 = 46732.0 / 5247,  a64 = 49.0 / 176, a65 = -5103.0 / 18656,
            b1  = 35.0 / 384,     b3  = 500.0 / 1113,    b4  = 125.0 / 192,     b5  = -2187.0 / 6784, b6 = 11.0 / 84,
            e1  = 71.0 / 57600,   e3  = -71.0 / 16695,   e4  = 71.0 / 1920,
            e5  = -17253.0 / 339200, e6 = 22.0 / 525,    e7 = -1.0 / 40;

        PacketState t, k2;
        for (int i = 0; i < NSTATE; ++i) t.y[i] = y.y[i] + h * (vsplat(a21) * k1.y[i]);
        geodesicRhsV(rc, t, k2);
        for (int i = 0; i < NSTATE; ++i) t.y[i] = y.y[i] + h * (vsplat(a31) * k1.y[i] + vsplat(a32) * k2.y[i]);
        geodesicRhsV(rc, t, k3);
        for (int i = 0; i < NSTATE; ++i) t.y[i] = y.y[i] + h * (vsplat(a41) * k1.y[i] + vsplat(a42) * k2.y[i] + vsplat(a43) * k3.y[i]);
        geodesicRhsV(rc, t, k4);
        for (int i = 0; i < NSTATE; ++i) t.y[i] = y.y[i] + h * (vsplat(a51) * k1.y[i] + vsplat(a52) * k2.y[i] + vsplat(a53) * k3.y[i] + vsplat(a54) * k4.y[i]);
        geodesicRhsV(rc, t, k5);
        for (int i = 0; i < NSTATE; ++i) t.y[i] = y.y[i] + h * (vsplat(a61) * k1.y[i] + vsplat(a62) * k2.y[i] + vsplat(a63) * k3.y[i] + vsplat(a64) * k4.y[i] + vsplat(a65) * k5.y[i]);
        geodesicRhsV(rc, t, k6);
        for (int i = 0; i < NSTATE; ++i)
            yOut.y[i] = y.y[i] + h * (vsplat(b1) * k1.y[i] + vsplat(b3) * k3.y[i] + vsplat(b4) * k4.y[i] + vsplat(b5) * k5.y[i] + vsplat(b6) * k6.y[i]);
        geodesicRhsV(rc, yOut, k7);

        vd err = vsplat(0.0);
        for (int i = 0; i < NSTATE; ++i) {
            vd e  = h * (vsplat(e1) * k1.y[i] + vsplat(e3) * k3.y[i] + vsplat(e4) * k4.y[i]
                       + vsplat(e5) * k5.y[i] + vsplat(e6) * k6.y[i] + vsplat(e7) * k7.y[i]);
            vd sc = atol + rtol * vmax(vabs(y.y[i]), vabs(yOut.y[i]));
            vd q  = e / sc;
            err = err + q * q;
        }
        return vsqrt(err * vsplat(1.0 / NSTATE));
    }

    void buildInterpolant(const PacketState& y0, const PacketState& y1, vd h) {
        constexpr double
            d1 = -12715105075.0 / 11282082432.0, d3 = 87487479700.0 / 32700410799.0,
            d4 = -10690763975.0 / 1880347072.0,  d5 = 701980252875.0 / 199316789632.0,
            d6 = -1453857185.0 / 822651844.0,    d7 = 69997945.0 / 29380423.0;
        for (int i = 0; i < NSTATE; ++i) {
            c1.y[i] = y0.y[i];
            c2.y[i] = y1.y[i] - y0.y[i];
            c3.y[i] = h * k1.y[i] - c2.y[i];
            c4.y[i] = c2.y[i] - h * k7.y[i] - c3.y[i];
            c5.y[i] = h * (vsplat(d1) * k1.y[i] + vsplat(d3) * k3.y[i] + vsplat(d4) * k4.y[i]
                         + vsplat(d5) * k5.y[i] + vsplat(d6) * k6.y[i] + vsplat(d7) * k7.y[i]);
        }
    }

    // Only theta is needed during the crossing search, so evaluate just that.
    vd interpTheta(vd s) const {
        vd s1 = vsplat(1.0) - s;
        return c1.y[1] + s * (c2.y[1] + s1 * (c3.y[1] + s * (c4.y[1] + s1 * c5.y[1])));
    }
    void interpolate(vd s, PacketState& out) const {
        vd s1 = vsplat(1.0) - s;
        for (int i = 0; i < NSTATE; ++i)
            out.y[i] = c1.y[i] + s * (c2.y[i] + s1 * (c3.y[i] + s * (c4.y[i] + s1 * c5.y[i])));
    }

    static vd nextStep(vd h, vd err, vd hMinL, vd hMaxL) {
        vd scale = vmin(vmax(vsplat(0.9) * vInvFifthRoot(err), vsplat(0.2)), vsplat(6.0));
        scale = vsel(err == vsplat(0.0), vsplat(6.0), scale);
        return vmin(vmax(h * scale, hMinL), hMaxL);
    }
};

// ---------------------------------------------------------------------------
// Packet propagation
// ---------------------------------------------------------------------------
struct PacketResult {
    Term term[LANES];
    State yEnd[LANES];
};

// Advance every active lane until it terminates.  Lanes that finish are parked
// on a benign state so their arithmetic stays finite and cannot slow the rest
// of the packet down; their results are already recorded.
inline void runPacket(const Kerr& kerr, const Disk& disk, const Propagator& prop,
                      const Geodesic g[LANES], const bool activeIn[LANES],
                      PacketResult& out, uint64_t* stepsOut = nullptr)
{
    PacketState y;
    vd E, L;
    vi alive = isplat(0);
    for (int i = 0; i < LANES; ++i) {
        bool on = activeIn[i];
        y.y[0][i] = on ? g[i].r   : 10.0;
        y.y[1][i] = on ? g[i].th  : HALF_PI;
        y.y[2][i] = on ? g[i].ph  : 0.0;
        y.y[3][i] = on ? g[i].pr  : 0.0;
        y.y[4][i] = on ? g[i].pth : 0.0;
        y.y[5][i] = 0.0;
        E[i]      = on ? g[i].E   : 1.0;
        L[i]      = on ? g[i].L   : 0.0;
        alive[i]  = on ? -1LL : 0LL;
        out.term[i] = Term::Captured;
    }
    if (!anyTrue(alive)) return;

    PacketConst rc = packConst(kerr.a, E, L);
    PacketDP dp;
    dp.rtol = vsplat(prop.rtol);

    geodesicRhsV(rc, y, dp.k1);

    const vd rEsc  = vsplat(prop.rEscape);
    const vd rCap  = vsplat(prop.rCapture);
    const vd rInV  = vsplat(disk.rIn);
    const vd rCutV = vsplat(disk.rCut);
    const vd zero  = vsplat(0.0);

    auto capFor = [&](const PacketState& s, const PacketState& k1) {
        vd c = vsplat(8.0);
        c = vmin(c, vsplat(0.30) * vmax(s.y[0], vsplat(1.0)) / vmax(vabs(k1.y[0]), vsplat(1e-30)));
        c = vmin(c, vsplat(0.25) / vmax(vabs(k1.y[1]), vsplat(1e-30)));
        c = vmin(c, vsplat(0.50) / vmax(vabs(k1.y[2]), vsplat(1e-30)));
        return c;
    };

    vd h = capFor(y, dp.k1);
    vi rejects = isplat(0);
    uint64_t steps = 0;

    // Retire a lane: record its outcome and park it.
    auto retire = [&](int i, Term t, const PacketState& src) {
        out.term[i] = t;
        for (int j = 0; j < NSTATE; ++j) out.yEnd[i].y[j] = src.y[j][i];
        alive[i] = 0;
        y.y[0][i] = 10.0; y.y[1][i] = HALF_PI; y.y[2][i] = 0.0;
        y.y[3][i] = 0.0;  y.y[4][i] = 0.0; y.y[5][i] = 0.0;
    };

    for (int step = 0; step < prop.maxSteps && anyTrue(alive); ++step) {
        ++steps;
        vd cap = capFor(y, dp.k1);
        h = vmin(h, cap);
        h = vmax(h, dp.hMin);

        PacketState yn;
        vd err = dp.trial(rc, y, h, yn);

        vi finite   = (err == err);                       // false only for NaN
        vi tooBig   = (err > vsplat(1.0)) | ~finite;
        vd hNew     = PacketDP::nextStep(h, err, dp.hMin, cap);
        vi stuck    = (rejects > isplat(40)) | (hNew <= dp.hMin * vsplat(1.001));
        vi accept   = alive & (~tooBig | stuck);
        vi reject   = alive & tooBig & ~stuck;

        rejects = vsel(reject, rejects + isplat(1), isplat(0));

        // Equatorial crossing, on accepted lanes only.
        vd t0 = vsplat(HALF_PI) - y.y[1];
        vd t1 = vsplat(HALF_PI) - yn.y[1];
        vi crossed = accept & ((t0 * t1) < zero);

        vi hitDisk = isplat(0);
        PacketState ym{};
        if (prop.diskActive && anyTrue(crossed)) {
            dp.buildInterpolant(y, yn, h);
            vd lo = zero, hi = vsplat(1.0);
            for (int it = 0; it < 40; ++it) {
                vd mid = vsplat(0.5) * (lo + hi);
                vd tm  = vsplat(HALF_PI) - dp.interpTheta(mid);
                vi keepLow = (tm * t0) > zero;
                lo = vsel(keepLow, mid, lo);
                hi = vsel(keepLow, hi, mid);
            }
            dp.interpolate(vsplat(0.5) * (lo + hi), ym);
            hitDisk = crossed & (ym.y[0] >= rInV) & (ym.y[0] <= rCutV);
        }

        // Advance accepted lanes that did not land on the disk.
        vi advance = accept & ~hitDisk;
        for (int i = 0; i < NSTATE; ++i) {
            y.y[i]     = vsel(advance, yn.y[i], y.y[i]);
            dp.k1.y[i] = vsel(advance, dp.k7.y[i], dp.k1.y[i]);
        }
        h = vsel(accept, hNew, vsel(reject, vsel(finite, hNew, h * vsplat(0.25)), h));

        // Retire lanes, scalar and rare.
        if (anyTrue(hitDisk)) {
            for (int i = 0; i < LANES; ++i)
                if (hitDisk[i] && alive[i]) retire(i, Term::Disk, ym);
        }
        vi escaped  = alive & (y.y[0] > rEsc);
        vi captured = alive & (y.y[3] < zero) & (y.y[0] < rCap);
        vi broken   = alive & ~(y.y[0] == y.y[0]);
        vi done = escaped | captured | broken;
        if (anyTrue(done)) {
            for (int i = 0; i < LANES; ++i) {
                if (!alive[i]) continue;
                if (escaped[i])       retire(i, Term::Escaped, y);
                else if (captured[i] || broken[i]) retire(i, Term::Captured, y);
            }
        }
    }

    // Anything still running hit the step budget.
    for (int i = 0; i < LANES; ++i)
        if (alive[i]) retire(i, Term::Stalled, y);

    if (stepsOut) *stepsOut += steps * uint64_t(LANES);
}

// ---------------------------------------------------------------------------
// Packet path tracing
// ---------------------------------------------------------------------------
//
// Propagation is vectorised; what happens at a surface is not.  Disk emission
// and scattering are branchy, per-lane and comparatively rare, so they stay
// scalar and reuse the same expressions as the scalar tracer.
inline void tracePacket(const Kerr& kerr, const Disk& disk, const Sky& sky,
                        const Propagator& prop, Geodesic g[LANES],
                        const Real lambda0[LANES], Rng rng[LANES],
                        int maxBounces, int nActive, Real radianceOut[LANES],
                        Real tObs = 0, uint64_t* steps = nullptr)
{
    Real travel[LANES] = {};      // light travel time back from the camera
    Real shift[LANES], throughput[LANES];
    int  bounce[LANES], crossings[LANES];
    bool alive[LANES];
    for (int i = 0; i < LANES; ++i) {
        radianceOut[i] = 0;
        shift[i] = 1; throughput[i] = 1;
        bounce[i] = 0; crossings[i] = 0;
        alive[i] = (i < nActive);
    }

    constexpr int MAX_CROSSINGS = 32;

    for (int segment = 0; segment < 64; ++segment) {
        bool any = false;
        for (int i = 0; i < LANES; ++i) any |= alive[i];
        if (!any) break;

        PacketResult res;
        runPacket(kerr, disk, prop, g, alive, res, steps);

        for (int i = 0; i < LANES; ++i) {
            if (!alive[i]) continue;
            Term t = res.term[i];

            if (t == Term::Captured || t == Term::Stalled) { alive[i] = false; continue; }

            if (t == Term::Escaped) {
                Real nuRatio = shift[i] * g[i].E;
                if (!(nuRatio > 0)) { alive[i] = false; continue; }
                Real gs = 1 / nuRatio;
                Real gs2 = gs * gs, gs5 = gs2 * gs2 * gs;
                Vec3 dir = asymptoticDirection(kerr, res.yEnd[i], g[i].E, g[i].L);
                radianceOut[i] += throughput[i] * gs5 *
                                  sky.radiance(dir, lambda0[i] * gs, bounce[i] == 0);
                alive[i] = false;
                continue;
            }

            // --- disk ---
            Real r   = res.yEnd[i].y[0];
            Real ph  = res.yEnd[i].y[2];
            Real pth = res.yEnd[i].y[4];

            travel[i] += res.yEnd[i].y[5];
            Real tEmit = tObs - travel[i];

            Vec4 u  = disk.orbitVelocity(r);
            Real om = disk.orbitOmega(r);
            Real nu = u[0] * (g[i].E - om * g[i].L);
            if (!(nu > 0)) { alive[i] = false; continue; }

            if (++crossings[i] > MAX_CROSSINGS) { alive[i] = false; continue; }

            Real mu = clampf(std::fabs(pth) / (r * nu), 1e-4, 1.0);
            Real trans = std::exp(-disk.opticalDepth(r) / mu);

            if (trans > 0 && rng[i].uniform() < trans) {
                g[i].r   = r;
                g[i].th  = HALF_PI + (pth > 0 ? 1e-11 : -1e-11);
                g[i].ph  = ph;
                g[i].pr  = res.yEnd[i].y[3];
                g[i].pth = pth;
                continue;                       // still alive, same momentum
            }

            shift[i] *= nu;
            Real gs = 1 / shift[i];
            Real gs2 = gs * gs, gs5 = gs2 * gs2 * gs;
            Real T = disk.temperature(r, ph, tEmit);
            if (T > 0) radianceOut[i] += throughput[i] * gs5 * spec::planck(lambda0[i] * gs, T);

            if (bounce[i]++ >= maxBounces) { alive[i] = false; continue; }
            if (rng[i].uniform() >= disk.albedo) { alive[i] = false; continue; }

            Real side = (pth > 0) ? -1.0 : 1.0;
            Metric m = kerr.metricAt(r, PI / 2);
            Vec4 e[4];
            buildTetrad(m, u, e);
            Vec3 s = cosineHemisphere(rng[i].uniform(), rng[i].uniform());
            Vec4 p;
            for (int j = 0; j < 4; ++j)
                p[j] = e[0][j] + s.x * e[1][j] + s.y * e[3][j] + (side * s.z) * e[2][j];
            Vec4 pl = m.lower(p);
            g[i].r   = r;
            g[i].th  = HALF_PI + side * 1e-11;
            g[i].ph  = ph;
            g[i].pr  = pl[1];
            g[i].pth = pl[2];
            g[i].E   = -pl[0];
            g[i].L   = pl[3];
            if (!(g[i].E > 0)) alive[i] = false;
        }
    }
}
