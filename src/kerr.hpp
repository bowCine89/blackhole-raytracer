// kerr.hpp -- Kerr spacetime: metric, null-geodesic flow, observer tetrads.
//
// Geometrised units with G = c = M = 1, so lengths are in gravitational radii
// and the spin parameter `a` is dimensionless in (-1, 1).
//
// Coordinates are Boyer-Lindquist (t, r, theta, phi).  The spacetime is
// stationary and axisymmetric, so E = -p_t and L = p_phi are constants of
// motion and only (r, theta, phi, p_r, p_theta) have to be integrated.
//
// We integrate in *Mino time* tau, defined by dtau = dlambda / Sigma.  With
//
//     F = Delta p_r^2 + p_theta^2 + (L - a E sin^2)^2 / sin^2
//                                 - [E (r^2 + a^2) - a L]^2 / Delta
//
// one has F = Sigma * g^{mu nu} p_mu p_nu, which vanishes on a null geodesic,
// and F/2 generates the flow in tau.  The Sigma factor cancels out of the r and
// theta equations entirely, which is why the right-hand side below is only a
// couple of dozen flops, branch-free, and free of any r-theta coupling.
#pragma once
#include "core.hpp"
#include <cstring>

struct Kerr {
    Real a = 0.0;

    explicit Kerr(Real spin = 0.0) : a(clampf(spin, -0.9999, 0.9999)) {}

    Real horizon() const { return 1 + std::sqrt(std::max<Real>(0, 1 - a * a)); }

    // Equatorial circular photon orbit.  A photon that falls inside the
    // prograde one while moving inward can never turn around again, which
    // gives a cheap and exact capture test that avoids integrating into the
    // stiff region where Delta -> 0.
    Real photonOrbit(bool prograde) const {
        Real s = prograde ? -a : a;
        return 2 * (1 + std::cos((2.0 / 3.0) * std::acos(clampf(s, -1, 1))));
    }

    // Innermost stable circular orbit (Bardeen, Press & Teukolsky 1972).
    Real isco(bool prograde) const {
        Real a2 = a * a;
        Real z1 = 1 + std::cbrt(1 - a2) * (std::cbrt(1 + a) + std::cbrt(1 - a));
        Real z2 = std::sqrt(3 * a2 + z1 * z1);
        Real t  = std::sqrt(std::max<Real>(0, (3 - z1) * (3 + z1 + 2 * z2)));
        return prograde ? 3 + z2 - t : 3 + z2 + t;
    }

    Real delta(Real r)           const { return r * r - 2 * r + a * a; }
    Real sigma(Real r, Real cth) const { return r * r + a * a * cth * cth; }
    Real bigA(Real r, Real sth2) const {              // (r^2+a^2)^2 - a^2 Delta sin^2
        Real rr = r * r + a * a;
        return rr * rr - a * a * delta(r) * sth2;
    }

    // Covariant metric at (r, theta).
    Metric metricAt(Real r, Real th) const {
        Real s = std::sin(th), c = std::cos(th);
        Real s2 = s * s;
        Real S = sigma(r, c), D = delta(r), A = bigA(r, s2);
        Metric m;
        m.g[0][0] = -(1 - 2 * r / S);
        m.g[0][3] = m.g[3][0] = -2 * a * r * s2 / S;
        m.g[1][1] = S / D;
        m.g[2][2] = S;
        m.g[3][3] = A * s2 / S;
        return m;
    }
};

// Phase-space point of a null geodesic.  E and L ride along as constants.
struct Geodesic {
    Real r = 0, th = 0, ph = 0;    // Boyer-Lindquist position (t is not needed)
    Real pr = 0, pth = 0;          // covariant momenta p_r, p_theta
    Real E = 1, L = 0;             // conserved -p_t and p_phi
};

// The five-component state the integrator advances.
// Six components: r, theta, phi, p_r, p_theta, t.
//
// Coordinate time rides along so that the *emission time* at a disk hit is
// known.  Light from the far side of the disk, and from the lensed images that
// loop around the hole, left earlier than light from the near side -- tens of M
// earlier, against an ISCO orbital period of ~24 M -- so without it an animated
// disk would show every part of itself at the same instant, which is wrong.
static constexpr int NSTATE = 6;
struct State { Real y[NSTATE]; };

// Everything about a ray that is fixed along its whole trajectory.  Hoisting
// these out of the right-hand side matters: the RHS runs six times per
// integration step and tens of times per pixel sample.
struct RayConst {
    Real a, a2, E, L;
    Real aL, aE, L2, a2E2;
    RayConst() = default;
    RayConst(Real spin, Real e, Real l)
        : a(spin), a2(spin * spin), E(e), L(l),
          aL(spin * l), aE(spin * e), L2(l * l), a2E2(spin * spin * e * e) {}
};

// Right-hand side of the Mino-time flow.  This is the hot loop of the whole
// renderer.  It is written to need exactly one sincos and two divisions: the
// naive form wants four divisions, and on this workload that is a measurable
// fraction of total render time.
inline void geodesicRhs(const RayConst& c, const State& s, State& out) {
    const Real r = s.y[0], th = s.y[1], pr = s.y[3], pth = s.y[4];

    Real sth, cth;
#if defined(__GNUC__) || defined(__clang__)
    __builtin_sincos(th, &sth, &cth);
#else
    sth = std::sin(th); cth = std::cos(th);
#endif
    // Keep the polar axis out of the denominators.  Geodesics with L != 0 have
    // a theta turning point before they reach it, so this clamp never fires on
    // a physical trajectory -- it only stops a stray trial step producing inf.
    if (std::fabs(sth) < 1e-8) sth = (sth < 0) ? -1e-8 : 1e-8;

    Real invS  = 1 / sth;
    Real invS2 = invS * invS;

    Real r2    = r * r;
    Real D     = r2 - 2 * r + c.a2;         // Delta
    Real dD    = 2 * r - 2;                 // dDelta/dr
    Real invD  = 1 / D;
    Real PinvD = (c.E * (r2 + c.a2) - c.aL) * invD;

    Real s2 = sth * sth;

    out.y[0] = D * pr;                                          // dr/dtau
    out.y[1] = pth;                                             // dtheta/dtau
    out.y[2] = c.L * invS2 - c.aE + c.a * PinvD;                // dphi/dtau
    // dp_r/dtau = -1/2 dF/dr
    out.y[3] = -0.5 * (dD * pr * pr - 4 * r * c.E * PinvD + PinvD * PinvD * dD);
    // dp_theta/dtau = -1/2 dF/dtheta
    out.y[4] = cth * (c.L2 * invS2 * invS - c.a2E2 * sth);
    // dt/dtau = -1/2 dF/dE, the same Hamiltonian, differentiated by energy.
    out.y[5] = c.a * (c.L - c.aE * s2) + (r2 + c.a2) * PinvD;
}

// Value of the null constraint F.  Zero analytically; its drift is a direct
// measure of integration error, which the --check mode reports.
inline Real nullConstraint(const Kerr& k, const State& s, Real E, Real L) {
    Real a = k.a, r = s.y[0], th = s.y[1];
    Real sth = std::sin(th); Real s2 = std::max<Real>(sth * sth, 1e-16);
    Real D = r * r - 2 * r + a * a;
    Real P = E * (r * r + a * a) - a * L;
    Real w = L - a * E * s2;
    return D * s.y[3] * s.y[3] + s.y[4] * s.y[4] + w * w / s2 - P * P / D;
}

// ---------------------------------------------------------------------------
// Dormand-Prince 5(4), adaptive, with the standard 5th-order dense output.
// The interpolant is what lets us locate the equatorial-plane crossing to
// near machine precision without a single extra force evaluation.
// ---------------------------------------------------------------------------
struct DormandPrince {
    Real rtol = 1e-9, atol = 1e-12;
    Real hMin = 1e-10, hMax = 8.0;

    // Stages kept from the last accepted step, for the dense output.
    State k1{}, k3{}, k4{}, k5{}, k6{}, k7{};
    State c1{}, c2{}, c3{}, c4{}, c5{};

    // One trial step of size h.  k1 must already hold f(y) (FSAL).
    // Returns the error norm; yOut and k7 are filled in either way.
    Real trial(const RayConst& rc, const State& y, Real h, State& yOut) {
        constexpr Real
            a21 = 1.0 / 5,
            a31 = 3.0 / 40,       a32 = 9.0 / 40,
            a41 = 44.0 / 45,      a42 = -56.0 / 15,      a43 = 32.0 / 9,
            a51 = 19372.0 / 6561, a52 = -25360.0 / 2187, a53 = 64448.0 / 6561,  a54 = -212.0 / 729,
            a61 = 9017.0 / 3168,  a62 = -355.0 / 33,     a63 = 46732.0 / 5247,  a64 = 49.0 / 176, a65 = -5103.0 / 18656,
            b1  = 35.0 / 384,     b3  = 500.0 / 1113,    b4  = 125.0 / 192,     b5  = -2187.0 / 6784, b6 = 11.0 / 84,
            e1  = 71.0 / 57600,   e3  = -71.0 / 16695,   e4  = 71.0 / 1920,
            e5  = -17253.0 / 339200, e6 = 22.0 / 525,    e7 = -1.0 / 40;

        State t, k2;
        for (int i = 0; i < NSTATE; ++i) t.y[i] = y.y[i] + h * a21 * k1.y[i];
        geodesicRhs(rc, t, k2);
        for (int i = 0; i < NSTATE; ++i) t.y[i] = y.y[i] + h * (a31 * k1.y[i] + a32 * k2.y[i]);
        geodesicRhs(rc, t, k3);
        for (int i = 0; i < NSTATE; ++i) t.y[i] = y.y[i] + h * (a41 * k1.y[i] + a42 * k2.y[i] + a43 * k3.y[i]);
        geodesicRhs(rc, t, k4);
        for (int i = 0; i < NSTATE; ++i) t.y[i] = y.y[i] + h * (a51 * k1.y[i] + a52 * k2.y[i] + a53 * k3.y[i] + a54 * k4.y[i]);
        geodesicRhs(rc, t, k5);
        for (int i = 0; i < NSTATE; ++i) t.y[i] = y.y[i] + h * (a61 * k1.y[i] + a62 * k2.y[i] + a63 * k3.y[i] + a64 * k4.y[i] + a65 * k5.y[i]);
        geodesicRhs(rc, t, k6);
        for (int i = 0; i < NSTATE; ++i)
            yOut.y[i] = y.y[i] + h * (b1 * k1.y[i] + b3 * k3.y[i] + b4 * k4.y[i] + b5 * k5.y[i] + b6 * k6.y[i]);
        geodesicRhs(rc, yOut, k7);      // FSAL: k1 of the next step

        Real err = 0;
        for (int i = 0; i < NSTATE; ++i) {
            Real e  = h * (e1 * k1.y[i] + e3 * k3.y[i] + e4 * k4.y[i] + e5 * k5.y[i] + e6 * k6.y[i] + e7 * k7.y[i]);
            Real sc = atol + rtol * std::max(std::fabs(y.y[i]), std::fabs(yOut.y[i]));
            Real q  = e / sc;
            err += q * q;
        }
        return std::sqrt(err / NSTATE);
    }

    // Dense-output coefficients for a step that was just accepted.
    void buildInterpolant(const State& y0, const State& y1, Real h) {
        constexpr Real
            d1 = -12715105075.0 / 11282082432.0, d3 = 87487479700.0 / 32700410799.0,
            d4 = -10690763975.0 / 1880347072.0,  d5 = 701980252875.0 / 199316789632.0,
            d6 = -1453857185.0 / 822651844.0,    d7 = 69997945.0 / 29380423.0;
        for (int i = 0; i < NSTATE; ++i) {
            c1.y[i] = y0.y[i];
            c2.y[i] = y1.y[i] - y0.y[i];
            c3.y[i] = h * k1.y[i] - c2.y[i];
            c4.y[i] = c2.y[i] - h * k7.y[i] - c3.y[i];
            c5.y[i] = h * (d1 * k1.y[i] + d3 * k3.y[i] + d4 * k4.y[i]
                         + d5 * k5.y[i] + d6 * k6.y[i] + d7 * k7.y[i]);
        }
    }

    // Evaluate the interpolant at s in [0,1] across the last accepted step.
    void interpolate(Real s, State& out) const {
        Real s1 = 1 - s;
        for (int i = 0; i < NSTATE; ++i)
            out.y[i] = c1.y[i] + s * (c2.y[i] + s1 * (c3.y[i] + s * (c4.y[i] + s1 * c5.y[i])));
    }

    // err^(-1/5) for the step-size controller.
    //
    // A general pow() call here measured at roughly 12% of total step time,
    // which is a lot for a quantity that only decides how big the next step is.
    // Seed from IEEE-754 exponent arithmetic, then three Newton iterations on
    // f(y) = y^-5 - x.  Measured worst case 2.6e-6 relative over the range the
    // controller uses, so the step sequence matches the exact version exactly.
    static Real invFifthRoot(Real x) {
        x = clampf(x, 1e-30, 1e30);
        uint64_t bits;
        std::memcpy(&bits, &x, 8);
        // bits(x^p) ~ p*(bits(x) - B) + B, with B tuned for the whole double.
        constexpr uint64_t B = 4607182418800017408ULL;   // 1023 << 52
        bits = uint64_t(int64_t(B) - (int64_t(bits) - int64_t(B)) / 5);
        Real y;
        std::memcpy(&y, &bits, 8);
        for (int i = 0; i < 3; ++i) {                    // y *= (6 - x*y^5)/5
            Real y2 = y * y, y5 = y2 * y2 * y;
            y *= (6.0 - x * y5) * 0.2;
        }
        return y;
    }

    static Real nextStep(Real h, Real err, Real hMinL, Real hMaxL) {
        Real scale = (err == 0) ? 6.0 : clampf(0.9 * invFifthRoot(err), 0.2, 6.0);
        return clampf(h * scale, hMinL, hMaxL);
    }
};

// ---------------------------------------------------------------------------
// Orthonormal tetrads
// ---------------------------------------------------------------------------

// Gram-Schmidt an orthonormal frame onto a normalised timelike u (g(u,u) = -1).
// Seeds are the coordinate directions r, theta, phi, so e[1], e[2], e[3] come
// out aligned with "radial", "polar" and "azimuthal" as closely as the metric
// allows.  e[0] is u itself.
inline void buildTetrad(const Metric& m, const Vec4& u, Vec4 e[4]) {
    e[0] = u;
    for (int n = 0; n < 3; ++n) {
        Vec4 v; v[n + 1] = 1;
        // Remove the timelike part first.  Because g(u,u) = -1, the projector
        // is v + g(v,u) u rather than v - g(v,u) u.
        Real pu = m.dot(v, e[0]);
        for (int i = 0; i < 4; ++i) v[i] += pu * e[0][i];
        for (int j = 1; j <= n; ++j) {
            Real pj = m.dot(v, e[j]);
            for (int i = 0; i < 4; ++i) v[i] -= pj * e[j][i];
        }
        Real nrm = std::sqrt(std::max<Real>(1e-300, m.dot(v, v)));
        for (int i = 0; i < 4; ++i) v[i] /= nrm;
        e[n + 1] = v;
    }
}

// Zero-angular-momentum observer (ZAMO) at (r, theta): the natural "locally
// non-rotating" camera frame, well defined all the way down through the
// ergosphere where static observers cease to exist.
struct Zamo {
    Real alpha;        // lapse
    Real omega;        // frame-dragging angular velocity
    Real sqrtSoverD;   // sqrt(Sigma/Delta)
    Real sqrtS;        // sqrt(Sigma)
    Real sSqrtAoverS;  // sin(theta) sqrt(A/Sigma)

    Zamo(const Kerr& k, Real r, Real th) {
        Real s = std::sin(th), c = std::cos(th), s2 = s * s;
        Real S = k.sigma(r, c), D = k.delta(r), A = k.bigA(r, s2);
        alpha       = std::sqrt(S * D / A);
        omega       = 2 * k.a * r / A;
        sqrtSoverD  = std::sqrt(S / D);
        sqrtS       = std::sqrt(S);
        sSqrtAoverS = s * std::sqrt(A / S);
    }

    // A photon of unit locally-measured energy travelling along the unit frame
    // direction n = (n_r, n_theta, n_phi) as measured by this observer.
    Geodesic emit(Real r, Real th, Real ph, const Vec3& n) const {
        Geodesic g;
        g.r = r; g.th = th; g.ph = ph;
        g.L   = n.z * sSqrtAoverS;
        g.E   = alpha + omega * g.L;
        g.pr  = n.x * sqrtSoverD;
        g.pth = n.y * sqrtS;
        return g;
    }
};
