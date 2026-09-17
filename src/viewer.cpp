// viewer.cpp -- interactive progressive viewer, SDL2 front end.
//
// The path tracer in render.hpp is the backend, unchanged: this file only
// decides *which* rays to trace and when to show the result.
//
// The whole design follows from one number.  At roughly 2.9 Mrays/s a 33 ms
// frame buys about 96k samples, while a 1280x720 window has 921k pixels.  Full
// resolution therefore cannot deliver even a tenth of a sample per pixel inside
// an interactive frame, so the viewer runs in two regimes:
//
//   moving  -- render at reduced resolution, one sample per pixel, each pixel
//              painted as a scale x scale block.  The scale is derived from
//              measured throughput to hit the frame budget.
//   settled -- full resolution, accumulating samples into a running mean
//              forever, so the image converges while you look at it.
//
// Workers never stop; they pull work items from a single monotonic counter and
// cycle through tiles indefinitely.  A pause barrier lets the main thread
// change camera or resolution safely: it flips `paused`, waits for the active
// worker count to reach zero, mutates, then releases.  Because no worker can be
// mid-tile at that moment, no stale sample is ever written into a fresh buffer
// and no per-tile generation tracking is needed.
#include "core.hpp"
#include "kerr.hpp"
#include "spectrum.hpp"
#include "scene.hpp"
#include "render.hpp"
#include "packet.hpp"
#include "image.hpp"
#include "video.hpp"

#include <SDL2/SDL.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include <memory>

using Clock = std::chrono::steady_clock;
static double since(Clock::time_point t) {
    return std::chrono::duration<double>(Clock::now() - t).count();
}
static inline Real fract(Real x) { return x - std::floor(x); }

// Tile size is chosen per resolution so there are always several tiles per
// worker; a fixed size would leave most threads idle at low render scales.
static int tileSizeFor(int rw, int rh, int threads) {
    double per = double(rw) * rh / std::max(1, 4 * threads);
    int t = int(std::sqrt(std::max(1.0, per)));
    return std::min(64, std::max(8, t));
}

// ---------------------------------------------------------------------------
// Scene and shared state
// ---------------------------------------------------------------------------
struct Scene {
    Kerr  kerr{0.94};
    Disk  disk;
    Sky   sky;
    Propagator prop;
    int   maxBounces = 2;
    uint64_t seed = 20260913;
    Real  lamScale = 1.0;      // LAMBDA_SPAN / integral of ybar
};

struct Target {
    int winW = 1280, winH = 720;   // window, and the size of `display`
    int scale = 4;                 // block size; 1 = full resolution
    int rw = 0, rh = 0;            // render resolution
};

struct Shared {
    std::mutex mtx;
    std::condition_variable cv;
    bool paused = false;
    // Set while the animation is playing back: there is nothing to trace, so the
    // workers park rather than burn every core on a finished image.
    bool idle   = false;
    bool quit   = false;
    int  activeWorkers = 0;

    // Mutated only while paused and activeWorkers == 0.
    Scene        scene;
    CameraParams cam;
    Target       rt;

    std::vector<double>   accum;    // rw*rh*3, CIE XYZ
    std::vector<uint32_t> count;    // rw*rh
    std::vector<uint8_t>  display;  // winW*winH*3, what workers paint into
    std::vector<uint8_t>  present;  // last *complete* pass, shown while moving

    int tilesX = 0, tilesY = 0, numTiles = 1, tileSize = 32;
    int nThreads = 1;

    // One flag per tile.  Work items w and w + numTiles name the same tile, so a
    // worker that falls a full pass behind could otherwise race another on the
    // same pixels.  Claiming the tile makes that impossible; a worker that finds
    // one taken simply drops that work item.
    std::unique_ptr<std::atomic<uint8_t>[]> tileBusy;

    // Tiles are visited in a fixed shuffled order rather than raster order.  A
    // partial pass then covers the frame uniformly instead of filling from the
    // top, which both makes the refinement read evenly and -- the reason it was
    // added -- lets the exposure meter trust a partly drawn frame.
    std::vector<int> tileOrder;

    std::atomic<uint64_t> work{0};          // monotonic work counter

    // While moving, dispatch is capped at exactly one item per tile.  Without
    // the cap the counter overshoots between main-thread checks, so some tiles
    // collect a second sample and appear visibly smoother than their
    // neighbours -- rectangular patches of differing noise.  0 means no cap,
    // which is what the settled regime wants.
    std::atomic<uint64_t> passLimit{0};
    std::atomic<uint64_t> samplesTraced{0};
    std::atomic<double>   exposure{1.0};
    std::atomic<uint32_t> minCount{0};      // samples per pixel, for the title
    double meterCovered = 0, meterAnchor = 0, meterTarget = 0;

    // Observer coordinate time, in M.  Mutated only under the pause barrier.
    Real tObs = 0.0;
};

// ---------------------------------------------------------------------------
// Tile rendering
// ---------------------------------------------------------------------------
static void renderTile(Shared& S, int tile, uint64_t sampleIdx) {
    const Scene&  sc = S.scene;
    const Target& rt = S.rt;

    Camera cam(sc.kerr, S.cam);

    int tx = (tile % S.tilesX) * S.tileSize;
    int ty = (tile / S.tilesX) * S.tileSize;
    int x1 = std::min(tx + S.tileSize, rt.rw);
    int y1 = std::min(ty + S.tileSize, rt.rh);

    const double expScale = S.exposure.load(std::memory_order_relaxed);
    uint64_t traced = 0;

    // Deposit one finished sample and repaint the block it stands for.
    auto shade = [&](size_t idx, int x, int y, const Real* lam, const Real* L) {
        double* a = &S.accum[idx * 3];
        for (int k = 0; k < spec::NLAMBDA; ++k) {
            if (L[k] <= 0) continue;
            Vec3 c = spec::cieXYZ(lam[k]);
            a[0] += c.x * L[k]; a[1] += c.y * L[k]; a[2] += c.z * L[k];
        }
        uint32_t n = ++S.count[idx];

        // Tone map this pixel now, while its data is hot in cache.
        Real inv = sc.lamScale / (Real(n) * spec::NLAMBDA) * expScale;
        Vec3 v = spec::xyzToLinearSrgb({a[0] * inv, a[1] * inv, a[2] * inv});
        v = spec::desaturateHighlights({std::max<Real>(0, v.x),
                                        std::max<Real>(0, v.y),
                                        std::max<Real>(0, v.z)}, 0.85);
        uint8_t r8 = uint8_t(clampf(spec::srgbEncode(spec::acesFilmic(v.x)) * 255.0 + 0.5, 0, 255));
        uint8_t g8 = uint8_t(clampf(spec::srgbEncode(spec::acesFilmic(v.y)) * 255.0 + 0.5, 0, 255));
        uint8_t b8 = uint8_t(clampf(spec::srgbEncode(spec::acesFilmic(v.z)) * 255.0 + 0.5, 0, 255));

        int bx0 = x * rt.scale, by0 = y * rt.scale;
        int bx1 = std::min(bx0 + rt.scale, rt.winW);
        int by1 = std::min(by0 + rt.scale, rt.winH);
        for (int by = by0; by < by1; ++by) {
            uint8_t* row = &S.display[(size_t(by) * rt.winW + bx0) * 3];
            for (int bx = bx0; bx < bx1; ++bx) { *row++ = r8; *row++ = g8; *row++ = b8; }
        }
    };

    // Eight horizontally adjacent pixels form a packet.  Neighbouring pixels at
    // the same sample index are the most coherent grouping available here, so
    // lane divergence stays low.
    for (int y = ty; y < y1; ++y) {
        for (int x = tx; x < x1; x += LANES) {
            const int n = std::min(LANES, x1 - x);

            Geodesic gp[LANES];
            Real     lam[LANES][spec::NLAMBDA];
            Rng      rgs[LANES];
            size_t   idxs[LANES];

            for (int j = 0; j < n; ++j) {
                const int px = x + j;
                const size_t idx = size_t(y) * rt.rw + px;
                idxs[j] = idx;

                // Low-discrepancy offsets: an R2 sequence in the sample index,
                // Cranley-Patterson rotated per pixel so neighbours decorrelate.
                uint32_t h = hashU32(uint32_t(idx) * 2654435761u + 1u);
                Real o1 = hashFloat(h);
                Real o2 = hashFloat(hashU32(h ^ 0x9e3779b9u));
                Real o3 = hashFloat(hashU32(h ^ 0x85ebca6bu));

                Real u1 = fract(0.7548776662 * Real(sampleIdx) + o1);
                Real u2 = fract(0.5698402910 * Real(sampleIdx) + o2);
                Real ul = fract(0.6180339887 * Real(sampleIdx) + o3);

                rgs[j] = Rng(idx + 1, sc.seed + sampleIdx * 0x9E3779B97F4A7C15ull);
                for (int k = 0; k < spec::NLAMBDA; ++k) {
                    Real f = ul + Real(k) / spec::NLAMBDA; f -= std::floor(f);
                    lam[j][k] = spec::LAMBDA_MIN + spec::LAMBDA_SPAN * f;
                }

                Real sx = 2.0 * (px + u1) / rt.rw - 1.0;
                Real sy = 1.0 - 2.0 * (y + u2) / rt.rh;
                gp[j] = cam.ray(sx, sy);
            }

            Real rad[LANES][spec::NLAMBDA];
            tracePacket(sc.kerr, sc.disk, sc.sky, sc.prop, gp, lam, rgs,
                        sc.maxBounces, n, rad, S.tObs);
            traced += uint64_t(n);

            for (int j = 0; j < n; ++j) shade(idxs[j], x + j, y, lam[j], rad[j]);
        }
    }
    S.samplesTraced.fetch_add(traced, std::memory_order_relaxed);
}

static void workerLoop(Shared& S) {
    for (;;) {
        {
            std::unique_lock<std::mutex> lk(S.mtx);
            S.cv.wait(lk, [&] { return (!S.paused && !S.idle) || S.quit; });
            if (S.quit) return;
            ++S.activeWorkers;
        }
        // Safe to touch scene/buffers unlocked: the main thread only mutates
        // them once activeWorkers has fallen to zero, which cannot happen
        // while we are counted.
        uint64_t limit = S.passLimit.load(std::memory_order_relaxed);
        if (limit && S.work.load(std::memory_order_relaxed) >= limit) {
            // Pass finished; wait for the main thread to publish and reset.
            { std::lock_guard<std::mutex> lk(S.mtx); --S.activeWorkers; }
            S.cv.notify_all();
            std::this_thread::sleep_for(std::chrono::microseconds(200));
            continue;
        }
        uint64_t w = S.work.fetch_add(1, std::memory_order_relaxed);
        if (limit && w >= limit) {
            { std::lock_guard<std::mutex> lk(S.mtx); --S.activeWorkers; }
            S.cv.notify_all();
            continue;
        }
        int      tile = S.tileOrder[size_t(w % uint64_t(S.numTiles))];
        uint64_t samp = w / uint64_t(S.numTiles);
        // Claim the tile.  If another worker is a full pass behind and still
        // inside this tile, drop the work item rather than race it; the only
        // consequence is that this tile skips one sample index, which the
        // per-pixel sample count already accounts for.
        if (S.tileBusy[tile].exchange(1, std::memory_order_acquire) == 0) {
            renderTile(S, tile, samp);
            S.tileBusy[tile].store(0, std::memory_order_release);
        }
        {
            std::lock_guard<std::mutex> lk(S.mtx);
            --S.activeWorkers;
        }
        S.cv.notify_all();
    }
}

// ---------------------------------------------------------------------------
// Reconfiguration behind the pause barrier
// ---------------------------------------------------------------------------
// `mutate` runs with every worker parked, so it may resize or clear anything.
template <typename F>
static void reconfigure(Shared& S, F&& mutate) {
    {
        std::unique_lock<std::mutex> lk(S.mtx);
        S.paused = true;
        S.cv.wait(lk, [&] { return S.activeWorkers == 0; });

        mutate();

        S.paused = false;
    }
    S.cv.notify_all();
}

// Resize the render buffers for the current scale and clear the accumulator.
// `display` is deliberately left alone: the previous image stays on screen and
// is overwritten tile by tile, which reads as a refinement rather than a flash
// of black.
static void resetAccumulation(Shared& S) {
    Target& rt = S.rt;
    rt.rw = std::max(1, (rt.winW + rt.scale - 1) / rt.scale);
    rt.rh = std::max(1, (rt.winH + rt.scale - 1) / rt.scale);
    S.tileSize = tileSizeFor(rt.rw, rt.rh, S.nThreads);
    S.tilesX = (rt.rw + S.tileSize - 1) / S.tileSize;
    S.tilesY = (rt.rh + S.tileSize - 1) / S.tileSize;
    S.numTiles = std::max(1, S.tilesX * S.tilesY);
    S.tileBusy = std::make_unique<std::atomic<uint8_t>[]>(size_t(S.numTiles));

    S.tileOrder.resize(size_t(S.numTiles));
    for (int i = 0; i < S.numTiles; ++i) S.tileOrder[i] = i;
    {   // deterministic shuffle, so runs stay reproducible
        Rng shuf(1234, 5678);
        for (int i = S.numTiles - 1; i > 0; --i) {
            int j = int(shuf.uniform() * (i + 1));
            std::swap(S.tileOrder[i], S.tileOrder[j <= i ? j : i]);
        }
    }

    S.accum.assign(size_t(rt.rw) * rt.rh * 3, 0.0);
    S.count.assign(size_t(rt.rw) * rt.rh, 0u);
    S.work.store(0, std::memory_order_relaxed);
}

// ---------------------------------------------------------------------------
// Exposure
// ---------------------------------------------------------------------------
// Sampled sparsely from the accumulator on the main thread.  Reading while
// workers write is a benign race: these are statistics, and an 8-byte aligned
// double load cannot tear on x86-64.  Smoothed over time so the picture does
// not pulse while samples arrive.
// Percentile of the luminance currently in the accumulator, plus how much of
// the frame that percentile was actually drawn from.
static double meterAccumulator(Shared& S, double& coverage) {
    const Target& rt = S.rt;
    static std::vector<double> lum;
    lum.clear();
    size_t n = size_t(rt.rw) * rt.rh;
    size_t stride = std::max<size_t>(1, n / 20000);
    size_t probed = 0;
    for (size_t i = 0; i < n; i += stride) {
        ++probed;
        if (S.count[i] == 0) continue;
        // Divide by NLAMBDA to match what the tone map actually displays; the
        // accumulator holds the sum over all wavelengths carried per path.
        lum.push_back(S.accum[i * 3 + 1] / (double(S.count[i]) * spec::NLAMBDA));
    }
    coverage = probed ? double(lum.size()) / double(probed) : 0.0;
    if (lum.empty()) return 0.0;
    // A slightly lower percentile than the batch renderer uses: at one sample
    // per pixel the top of the distribution is dominated by outliers.
    size_t k = std::min(lum.size() - 1, size_t(lum.size() * 0.99));
    std::nth_element(lum.begin(), lum.begin() + k, lum.end());
    return lum[k];
}

// Exposure is a *ratio*, so it is carried and smoothed in log space.  The old
// linear form let one bad reading throw it to 1e15 and then took dozens of
// frames to crawl back, which is what made the first pass after a move glare.
static double updateExposure(Shared& S, double key, double userExposure, double logSmoothed) {
    double coverage = 0;
    double anchor = meterAccumulator(S, coverage);
    S.meterCovered = coverage;
    S.meterAnchor = anchor;

    // Never meter from a frame that is barely drawn.  Partial passes used to
    // produce anchors 5 orders of magnitude too small; holding the last good
    // reading is strictly better than believing a bad one.
    if (coverage < 0.25 || !(anchor > 0)) {
        S.exposure.store(std::exp(logSmoothed) * userExposure, std::memory_order_relaxed);
        return logSmoothed;
    }

    double target = key / (anchor * S.scene.lamScale);
    target = clampf(target, 1e-6, 1e6);
    S.meterTarget = target;

    double logTarget = std::log(target);
    if (!(logSmoothed > -1e30 && logSmoothed < 1e30)) logSmoothed = logTarget;   // first reading
    logSmoothed += (logTarget - logSmoothed) * 0.25;
    S.exposure.store(std::exp(logSmoothed) * userExposure, std::memory_order_relaxed);
    return logSmoothed;
}

// One coarse synchronous pass before anything is displayed, so frame one is
// already correctly exposed instead of arriving white and settling down.
static double meterInitial(const Scene& sc, const CameraParams& cam, double key) {
    Camera c(sc.kerr, cam);
    const int w = 96, h = 54, spp = 4;
    std::vector<double> lum;
    lum.reserve(size_t(w) * h);
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            double acc = 0;
            for (int s = 0; s < spp; ++s) {
                uint32_t hsh = hashCombine(uint32_t(x), uint32_t(y), uint32_t(s) + 1u);
                Real u1 = hashFloat(hsh), u2 = hashFloat(hashU32(hsh ^ 0x9e3779b9u));
                Real ul = fract(0.6180339887 * s + hashFloat(hashU32(hsh ^ 0x85ebca6bu)));
                Rng rng(size_t(y) * w + x + 1, sc.seed + uint64_t(s) * 0x9E3779B97F4A7C15ull);
                Real lam[spec::NLAMBDA];
                for (int k = 0; k < spec::NLAMBDA; ++k) {
                    Real f = ul + Real(k) / spec::NLAMBDA; f -= std::floor(f);
                    lam[k] = spec::LAMBDA_MIN + spec::LAMBDA_SPAN * f;
                }
                Real sx = 2.0 * (x + u1) / w - 1.0;
                Real sy = 1.0 - 2.0 * (y + u2) / h;
                Geodesic g = c.ray(sx, sy);
                Real rad[spec::NLAMBDA];
                tracePath(sc.kerr, sc.disk, sc.sky, sc.prop, g, lam, rad, rng, sc.maxBounces);
                for (int k = 0; k < spec::NLAMBDA; ++k)
                    if (rad[k] > 0) acc += spec::cieXYZ(lam[k]).y * rad[k];
            }
            lum.push_back(acc / (spp * spec::NLAMBDA));
        }
    }
    size_t k = std::min(lum.size() - 1, size_t(lum.size() * 0.99));
    std::nth_element(lum.begin(), lum.begin() + k, lum.end());
    double anchor = lum[k];
    if (!(anchor > 0)) return std::log(1.0);
    return std::log(clampf(key / (anchor * sc.lamScale), 1e-6, 1e6));
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Bake-and-loop
// ---------------------------------------------------------------------------
//
// The viewer cannot animate live -- a turning disk changes faster than the
// tracer converges -- but it can *bake*: hold the camera still, converge one
// frame, advance the clock, repeat, then loop what it captured.
//
// The interface is one context-sensitive key.  Space always does the obvious
// next thing, and the window title always says what that is:
//
//     LIVE     -> space starts baking
//     BAKING   -> space stops and plays what has been captured
//     PLAYING  -> space returns to live
//
// Escape backs out of whatever is happening without quitting.
enum class Mode { Live, Baking, Playing, Rays, Hero };

static const char* modeName(Mode m) {
    switch (m) {
        case Mode::Baking:  return "BAKING";
        case Mode::Playing: return "PLAYING";
        case Mode::Rays:    return "RAYS";
        case Mode::Hero:    return "HERO";
        default:            return "live";
    }
}

// ---------------------------------------------------------------------------
// Ray visualiser
// ---------------------------------------------------------------------------
//
// G stops the tracer and shows the geometry it was working in.  A sparse grid
// of the camera's own rays is integrated once and kept as polylines; dragging
// then moves *your* eye around them, not theirs -- the rays stay pinned to the
// viewpoint they were launched from, which is the whole point.  Watching a ray
// wind twice around the photon ring and come back out explains the picture in a
// way the picture cannot.
//
// The rays are traced backwards, as the renderer traces them: they leave the
// camera and run out into the spacetime.  Colour is the outcome, which is
// exactly the classification that makes the final image what it is.

// Boyer-Lindquist to Cartesian, in the oblate-spheroidal embedding that goes
// with these coordinates: surfaces of constant r are ellipsoids of semi-axis
// sqrt(r^2 + a^2) across and r along the spin axis, so the horizon shows its
// true flattening rather than being drawn as a sphere it is not.
static Vec3 blToCart(Real r, Real th, Real ph, Real a) {
    Real s = std::sin(th), c = std::cos(th);
    Real rho = std::sqrt(r * r + a * a);
    return { rho * s * std::cos(ph), rho * s * std::sin(ph), r * c };
}

struct RayLine {
    std::vector<Vec3> pts;
    Term end = Term::Escaped;
};

struct RayVis {
    std::vector<RayLine> lines;
    Vec3   target{0, 0, 0};       // what the inspection camera orbits
    CameraParams source;          // the viewpoint the rays were launched from
    Real   sourceA = 0;
    Real   rHorizon = 2, rIn = 6, rOut = 18;
    int    gridX = 21, gridY = 13;
    bool   captured = false;

    // The traced frame, frozen on entry and hung in space where the camera's
    // image plane actually is.  Seeing the picture the rays produced, sitting
    // in the bundle that produced it, is the whole reason this is worth having.
    std::vector<uint8_t> imgPixels;      // winW*winH*3 at the moment of capture
    int          imgW = 0, imgH = 0;
    SDL_Texture* imgTex = nullptr;
    Vec3 camPos{}, camFwd{}, camRight{}, camUp{};
    Real planeDist = 0, planeHalfW = 0, planeHalfH = 0;
};

// Cartesian images of the Boyer-Lindquist basis directions, by finite
// difference on the embedding.  Cheap, and it keeps the plane's orientation
// tied to the same basis the camera is built from rather than to a guess.
static void cartBasis(Real r, Real th, Real ph, Real a,
                      Vec3& er, Vec3& eth, Vec3& eph) {
    const Real h = 1e-4;
    Vec3 p = blToCart(r, th, ph, a);
    er  = normalize(blToCart(r + h, th, ph, a) - p);
    eth = normalize(blToCart(r, th + h, ph, a) - p);
    eph = normalize(blToCart(r, th, ph + h, a) - p);
}

// An ordinary pinhole camera in the embedding space.  Nothing here is a
// geodesic: we are looking *at* the spacetime from outside it, so straight
// lines are the honest choice and curvature belongs to what is being drawn.
struct ViewCam {
    Vec3 eye, fwd, right, up;
    Real tanHalf = 1, aspect = 1;
    ViewCam(const CameraParams& p, Real aspectIn, Vec3 target = {0, 0, 0}) {
        Real th = p.incDeg * PI / 180, ph = p.phiDeg * PI / 180;
        eye = target + Vec3{ p.camR * std::sin(th) * std::cos(ph),
                             p.camR * std::sin(th) * std::sin(ph),
                             p.camR * std::cos(th) };
        fwd = normalize(target - eye);
        Vec3 wup{0, 0, 1};
        if (std::fabs(dot(fwd, wup)) > 0.999) wup = Vec3{0, 1, 0};
        right = normalize(cross(fwd, wup));
        up    = normalize(cross(right, fwd));
        tanHalf = std::tan(0.5 * clampf(p.fovDeg, 5.0, 120.0) * PI / 180);
        aspect  = aspectIn;
    }
    // Returns false behind the eye; depth comes back in z for near-plane clipping.
    bool project(const Vec3& p, Real& sx, Real& sy, Real& z) const {
        Vec3 d = p - eye;
        z = dot(d, fwd);
        if (z <= 1e-4) return false;
        sx = dot(d, right) / (z * tanHalf * aspect);
        sy = dot(d, up)    / (z * tanHalf);
        return true;
    }
};

// One segment, clipped against the near plane so a polyline that passes behind
// the eye does not wrap across the screen.
static void drawSeg(SDL_Renderer* ren, const ViewCam& vc, int w, int h,
                    Vec3 a, Vec3 b) {
    const Real near = 1e-3;
    Real za = dot(a - vc.eye, vc.fwd), zb = dot(b - vc.eye, vc.fwd);
    if (za <= near && zb <= near) return;
    if (za <= near) { Real t = (near - za) / (zb - za); a = a + (b - a) * t; }
    if (zb <= near) { Real t = (near - zb) / (za - zb); b = b + (a - b) * t; }
    Real ax, ay, az, bx, by, bz;
    if (!vc.project(a, ax, ay, az) || !vc.project(b, bx, by, bz)) return;
    auto px = [&](Real sx) { return int((sx + 1) * 0.5 * w); };
    auto py = [&](Real sy) { return int((1 - sy) * 0.5 * h); };
    SDL_RenderDrawLine(ren, px(ax), py(ay), px(bx), py(by));
}

static void drawCircle(SDL_Renderer* ren, const ViewCam& vc, int w, int h,
                       Real radius, Real zHeight, int segs = 96) {
    Vec3 prev{radius, 0, zHeight};
    for (int i = 1; i <= segs; ++i) {
        Real t = TWO_PI * i / segs;
        Vec3 cur{radius * std::cos(t), radius * std::sin(t), zHeight};
        drawSeg(ren, vc, w, h, prev, cur);
        prev = cur;
    }
}

// Trace the grid once.  Escapes are cut off a little beyond the camera so the
// picture stays about the hole rather than about a sphere of radius 2000.
static void captureRays(RayVis& rv, const Scene& sc, const CameraParams& cam) {
    rv.lines.clear();
    rv.source   = cam;
    rv.sourceA  = sc.kerr.a;
    rv.rHorizon = sc.kerr.horizon();
    rv.rIn      = sc.disk.rIn;
    rv.rOut     = sc.disk.rOut;
    rv.target   = Vec3{0, 0, 0};    // orbit the hole; entry sits at the launch point

    // The image plane, in the same frame the renderer's Camera uses: fwd is
    // -e_r, up is -e_theta, right is +e_phi, then yaw and pitch on top.  Its
    // half-extents follow tanHalf exactly, so the quad is the frustum the
    // frame was actually rendered through and not an approximation of it.
    Vec3 er, eth, eph;
    cartBasis(cam.camR, cam.incDeg * PI / 180, cam.phiDeg * PI / 180, sc.kerr.a, er, eth, eph);
    rv.camPos   = blToCart(cam.camR, cam.incDeg * PI / 180, cam.phiDeg * PI / 180, sc.kerr.a);
    rv.camFwd   = -er;
    rv.camUp    = -eth;
    rv.camRight = eph;
    Camera::rotate(rv.camUp,    cam.yawDeg   * PI / 180, rv.camFwd, rv.camRight);
    Camera::rotate(rv.camRight, cam.pitchDeg * PI / 180, rv.camFwd, rv.camUp);
    Real tanHalf   = std::tan(0.5 * cam.fovDeg * PI / 180);
    rv.planeDist   = cam.camR * 0.42;
    rv.planeHalfH  = rv.planeDist * tanHalf;
    rv.planeHalfW  = rv.planeHalfH * cam.aspect;

    Camera camera(sc.kerr, cam);
    Propagator prop = sc.prop;
    prop.rEscape = cam.camR * 2.4;
    prop.maxSteps = 20000;

    std::vector<Vec3> bl;
    for (int j = 0; j < rv.gridY; ++j)
        for (int i = 0; i < rv.gridX; ++i) {
            Real sx = 2.0 * (i + 0.5) / rv.gridX - 1.0;
            Real sy = 1.0 - 2.0 * (j + 0.5) / rv.gridY;
            Geodesic g = camera.ray(sx, sy);
            State yEnd{};
            bl.clear();
            Term t = prop.run(g, yEnd, nullptr, &bl);

            RayLine rl;
            rl.end = t;
            rl.pts.reserve(bl.size() + 1);
            for (const Vec3& q : bl) rl.pts.push_back(blToCart(q.x, q.y, q.z, rv.sourceA));
            rl.pts.push_back(blToCart(yEnd.y[0], yEnd.y[1], yEnd.y[2], rv.sourceA));
            rv.lines.push_back(std::move(rl));
        }
    rv.captured = true;
}

// Interpolate two inspection poses.  Azimuth takes the short way round, so a
// return that crosses phi = 0 does not unwind the long way.
static CameraParams lerpPose(const CameraParams& a, const CameraParams& b, Real u) {
    CameraParams o = a;
    Real dPhi = b.phiDeg - a.phiDeg;
    while (dPhi >  180) dPhi -= 360;
    while (dPhi < -180) dPhi += 360;
    o.phiDeg   = a.phiDeg + dPhi * u;
    o.incDeg   = a.incDeg + (b.incDeg - a.incDeg) * u;
    o.fovDeg   = a.fovDeg + (b.fovDeg - a.fovDeg) * u;
    o.yawDeg   = a.yawDeg + (b.yawDeg - a.yawDeg) * u;
    o.pitchDeg = a.pitchDeg + (b.pitchDeg - a.pitchDeg) * u;
    // Radius moves geometrically: a linear sweep from 300M to 40M spends most
    // of its time far away and then lunges at the end.
    o.camR = a.camR * std::pow(b.camR / std::max<Real>(a.camR, 1e-6), double(u));
    return o;
}

static void rebuildRayTexture(SDL_Renderer* ren, RayVis& rv) {
    if (rv.imgTex) { SDL_DestroyTexture(rv.imgTex); rv.imgTex = nullptr; }
    if (rv.imgW <= 0 || rv.imgH <= 0 ||
        rv.imgPixels.size() < size_t(rv.imgW) * size_t(rv.imgH) * 3) return;
    rv.imgTex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_RGB24,
                                  SDL_TEXTUREACCESS_STATIC, rv.imgW, rv.imgH);
    if (!rv.imgTex) return;
    SDL_UpdateTexture(rv.imgTex, nullptr, rv.imgPixels.data(), rv.imgW * 3);
    SDL_SetTextureBlendMode(rv.imgTex, SDL_BLENDMODE_NONE);
}

static void drawRayVis(SDL_Renderer* ren, const RayVis& rv,
                       const CameraParams& view, int w, int h) {
    SDL_SetRenderDrawColor(ren, 6, 7, 10, 255);
    SDL_RenderClear(ren);
    ViewCam vc(view, Real(w) / Real(std::max(1, h)), rv.target);

    // Equatorial grid: the disk between rIn and rOut, then a few rings beyond
    // it for scale.  Drawn first so the rays read on top of it.
    SDL_SetRenderDrawColor(ren, 122, 88, 44, 255);
    for (int i = 0; i <= 6; ++i)
        drawCircle(ren, vc, w, h, rv.rIn + (rv.rOut - rv.rIn) * i / 6.0, 0.0);
    SDL_SetRenderDrawColor(ren, 52, 52, 64, 255);
    for (int k = 1; k <= 3; ++k)
        drawCircle(ren, vc, w, h, rv.rOut + k * (rv.rOut - rv.rIn) * 0.5, 0.0);
    // Spokes, so rotation about the axis is readable.
    SDL_SetRenderDrawColor(ren, 92, 68, 36, 255);
    for (int k = 0; k < 12; ++k) {
        Real t = TWO_PI * k / 12;
        drawSeg(ren, vc, w, h,
                Vec3{rv.rIn * std::cos(t), rv.rIn * std::sin(t), 0},
                Vec3{rv.rOut * std::cos(t), rv.rOut * std::sin(t), 0});
    }

    // The horizon, as the oblate surface it actually is.
    SDL_SetRenderDrawColor(ren, 120, 120, 140, 255);
    const Real rh = rv.rHorizon, a = rv.sourceA;
    for (int m = 0; m < 12; ++m) {                       // meridians
        Real ph = PI * m / 12;
        Vec3 prev = blToCart(rh, 0.0, ph, a);
        for (int i = 1; i <= 48; ++i) {
            Vec3 cur = blToCart(rh, PI * i / 48, ph, a);
            drawSeg(ren, vc, w, h, prev, cur);
            prev = cur;
        }
    }
    for (int l = 1; l < 8; ++l) {                        // parallels
        Real th = PI * l / 8;
        Real rr = std::sqrt(rh * rh + a * a) * std::sin(th);
        drawCircle(ren, vc, w, h, rr, rh * std::cos(th), 48);
    }

    // The traced frame, hung at the image plane.  Drawn before the rays so the
    // bundle reads on top of it; with no depth buffer that is the better of the
    // two wrong answers, and most of a ray's length lies beyond the plane anyway.
    if (rv.imgTex && rv.imgW > 0 && rv.imgH > 0) {
        // RenderGeometry maps affinely per triangle, so an oblique plane needs
        // subdividing or the seams show.  This is dense enough that they do not.
        const int NU = 40, NV = 26;
        auto corner = [&](int i, int j) {
            Real u = Real(i) / NU, v = Real(j) / NV;
            return rv.camPos + rv.camFwd * rv.planeDist
                 + rv.camRight * ((2 * u - 1) * rv.planeHalfW)
                 + rv.camUp    * ((1 - 2 * v) * rv.planeHalfH);
        };
        auto toPix = [&](const Vec3& q, SDL_FPoint& out) {
            Real sx, sy, z;
            if (!vc.project(q, sx, sy, z)) return false;
            out.x = float((sx + 1) * 0.5 * w);
            out.y = float((1 - sy) * 0.5 * h);
            return true;
        };
        std::vector<SDL_Vertex> verts;
        std::vector<int> idx;
        verts.reserve(size_t(NU + 1) * (NV + 1));
        std::vector<char> ok(size_t(NU + 1) * (NV + 1), 0);
        const SDL_Color tint{255, 255, 255, 255};
        for (int j = 0; j <= NV; ++j)
            for (int i = 0; i <= NU; ++i) {
                SDL_Vertex v{};
                v.color = tint;
                v.tex_coord = { float(i) / NU, float(j) / NV };
                ok[size_t(j) * (NU + 1) + i] = toPix(corner(i, j), v.position) ? 1 : 0;
                verts.push_back(v);
            }
        auto at = [&](int i, int j) { return size_t(j) * (NU + 1) + i; };
        for (int j = 0; j < NV; ++j)
            for (int i = 0; i < NU; ++i) {
                size_t a0 = at(i, j), b0 = at(i + 1, j), c0 = at(i, j + 1), d0 = at(i + 1, j + 1);
                if (!(ok[a0] && ok[b0] && ok[c0] && ok[d0])) continue;   // straddles the eye
                idx.push_back(int(a0)); idx.push_back(int(b0)); idx.push_back(int(d0));
                idx.push_back(int(a0)); idx.push_back(int(d0)); idx.push_back(int(c0));
            }
        if (!idx.empty())
            SDL_RenderGeometry(ren, rv.imgTex, verts.data(), int(verts.size()),
                               idx.data(), int(idx.size()));
        // A border, so the plane reads as an object rather than a smear.
        SDL_SetRenderDrawColor(ren, 150, 150, 165, 255);
        Vec3 c00 = corner(0, 0), c10 = corner(NU, 0), c11 = corner(NU, NV), c01 = corner(0, NV);
        drawSeg(ren, vc, w, h, c00, c10); drawSeg(ren, vc, w, h, c10, c11);
        drawSeg(ren, vc, w, h, c11, c01); drawSeg(ren, vc, w, h, c01, c00);
        // And the four frustum edges back to the eye.
        SDL_SetRenderDrawColor(ren, 70, 70, 86, 255);
        drawSeg(ren, vc, w, h, rv.camPos, c00); drawSeg(ren, vc, w, h, rv.camPos, c10);
        drawSeg(ren, vc, w, h, rv.camPos, c11); drawSeg(ren, vc, w, h, rv.camPos, c01);
    }

    // The rays, coloured by how each one ended.
    for (const RayLine& rl : rv.lines) {
        switch (rl.end) {
            case Term::Captured: SDL_SetRenderDrawColor(ren, 232,  72,  60, 255); break;
            case Term::Disk:     SDL_SetRenderDrawColor(ren, 255, 176,  64, 255); break;
            case Term::Escaped:  SDL_SetRenderDrawColor(ren, 120, 190, 255, 255); break;
            default:             SDL_SetRenderDrawColor(ren, 110, 110, 110, 255); break;
        }
        for (size_t i = 1; i < rl.pts.size(); ++i)
            drawSeg(ren, vc, w, h, rl.pts[i - 1], rl.pts[i]);
    }

    // Where the rays came from.
    SDL_SetRenderDrawColor(ren, 255, 255, 255, 255);
    Vec3 src = blToCart(rv.source.camR, rv.source.incDeg * PI / 180,
                        rv.source.phiDeg * PI / 180, rv.sourceA);
    Real mx, my, mz;
    if (vc.project(src, mx, my, mz)) {
        int cx = int((mx + 1) * 0.5 * w), cy = int((1 - my) * 0.5 * h);
        for (int d = -5; d <= 5; ++d) {
            SDL_RenderDrawPoint(ren, cx + d, cy);
            SDL_RenderDrawPoint(ren, cx, cy + d);
        }
    }
}

// ---------------------------------------------------------------------------
// The hero shot
// ---------------------------------------------------------------------------
//
// A scripted move, rendered at whatever resolution is asked for rather than at
// the window's, and streamed to the video a frame at a time.  Streaming is not
// an optimisation here: a minute of 4K is 1440 frames of 24 MB, and the bake's
// habit of keeping every frame in memory would want 35 GB of it.
//
// The move starts below the disk plane looking up, framed on the approaching
// side -- which is the bright one, by a factor of about four at the default
// inclination -- then pushes in and rises through the plane, settling centred
// on the shadow.  Crossing happens around r = 22 M, comfortably outside the
// disk's outer edge at 18 M, so the camera passes over the rim rather than
// through it.
struct Hero {
    bool   active = false;
    int    w = 3840, h = 2160;      // 4K by default
    int    frames = 1440;           // 60 s at 24 fps
    int    idx = 0;
    int    spp = 96;                // convergence before the clock advances
    double fps = 24.0;
    Real   tSpan = 48.0;            // coordinate time covered, about two ISCO orbits
    Real   t0 = 0;
    std::string path;
    vid::VideoWriter vw;
    bool   videoOpen = false;
    double lockedExposure = 0;

    // what to put back afterwards
    int          saveW = 0, saveH = 0, saveScale = 1;
    CameraParams saveCam;
    Real         saveT = 0;
    Clock::time_point began{};
};

// Position along the move, u running 0 to 1.  Smoothstep so it eases in and
// out rather than starting and stopping abruptly; the radius moves
// geometrically, because a linear sweep from 46 M to 11 M spends most of its
// time far away and then lunges at the end.
static CameraParams heroPose(const CameraParams& base, double u) {
    double e = u * u * (3 - 2 * u);
    CameraParams p = base;
    p.incDeg   = 106.0 + (72.0 - 106.0) * e;      // below the plane, then above it
    // 46 M to 24 M, not further.  The shadow is about 5 M across and at 11 M it
    // overflows the frame entirely; ending at 24 M puts it at roughly two
    // thirds of the half-height, which fills the frame without losing its shape.
    p.camR     = Real(46.0 * std::pow(24.0 / 46.0, e));
    p.fovDeg   = 42.0 + (34.0 - 42.0) * e;
    p.yawDeg   = Real(-7.0 * (1.0 - e));          // negative aims right, at the glare
    p.pitchDeg = 0;
    return p;
}

struct Bake {
    std::vector<std::vector<uint8_t>> frames;   // tone-mapped, window-sized
    int    targetSpp = 64;       // convergence before the clock advances
    Real   step      = -1;       // M between frames; <0 = one ISCO orbit / 96
    int    maxFrames = 180;      // ~500 MB at 720p
    Real   t0        = 0;        // clock at the first frame

    double fps       = 24.0;     // playback rate
    double playHead  = 0;        // fractional frame index
    double lockedExposure = 0;   // frozen during a bake, so playback cannot flicker

    Real   orbitPeriod = 24.0;   // ISCO orbital period, the natural loop unit
    Real   loopSpan    = 0;      // planned sequence length; the loop closes on it
    std::string videoPath;       // empty disables the video
    int         crf = 18;        // H.265 quality when videoPath is .mp4/.mkv

    size_t bytes() const { return frames.empty() ? 0 : frames.size() * frames[0].size(); }
    Real   span()  const { return frames.empty() ? 0 : Real(frames.size()) * step; }
    Real   orbits() const { return span() / orbitPeriod; }
};

// Write the captured sequence as an uncompressed AVI: lossless, no dependency,
// and directly playable.  Each frame cost seconds to converge, so there is no
// case for throwing quality away on the way to disk.
static bool writeBakeVideo(const Bake& bake, int w, int h) {
    if (bake.videoPath.empty() || bake.frames.size() < 2) return false;
    vid::VideoWriter vw;
    vw.crf = bake.crf;
    if (!vw.open(bake.videoPath, w, h, bake.fps)) {
        std::printf("  could not open %s for writing\n", bake.videoPath.c_str());
        return false;
    }
    for (const auto& f : bake.frames) vw.addFrame(f.data());
    if (!vw.close())
        std::printf("  the encoder reported a failure; %s may be truncated\n",
                    bake.videoPath.c_str());
    double mb = vw.bytes() / 1.0e6;
    char size[64];
    if (mb >= 1.0) std::snprintf(size, sizeof size, "%.1f MB", mb);
    else           std::snprintf(size, sizeof size, "%.0f KB", vw.bytes() / 1024.0);
    std::printf("  wrote %s   %d frames, %dx%d, %.0f fps, %s (%s)   (%.2f ISCO orbits)\n",
                bake.videoPath.c_str(), vw.frames(), w, h, bake.fps,
                size, vw.codec(), double(bake.orbits()));
    std::fflush(stdout);
    return true;
}

// Read the rendered view straight back off the renderer.  The ray visualiser
// draws with SDL primitives rather than into a pixel buffer, so this is the only
// way to get a picture of it out.
static bool savePresented(SDL_Renderer* ren, int w, int h, const std::string& path) {
    std::vector<uint8_t> rgba(size_t(w) * h * 4);
    if (SDL_RenderReadPixels(ren, nullptr, SDL_PIXELFORMAT_ABGR8888,
                             rgba.data(), w * 4) != 0)
        return false;
    std::vector<uint8_t> rgb(size_t(w) * h * 3);
    for (size_t i = 0; i < size_t(w) * h; ++i) {
        rgb[i * 3 + 0] = rgba[i * 4 + 0];
        rgb[i * 3 + 1] = rgba[i * 4 + 1];
        rgb[i * 3 + 2] = rgba[i * 4 + 2];
    }
    return img::writePng(path, w, h, rgb);
}

static void printControls() {
    std::printf(
"\nControls\n"
"  SPACE                bake an animation -> play it -> back to live\n"
"  ESC                  back out of baking or playback\n"
"\n"
"  left drag            orbit (inclination and azimuth)\n"
"  right drag           pan the aim point\n"
"  wheel                dolly in/out\n"
"  ctrl + wheel         field of view\n"
"  [ / ]                spin a/M down / up\n"
"  , / .                disk outer radius   (playback speed while playing)\n"
"  - / =                exposure\n"
"  b                    cycle scattering bounces (0-4)\n"
"  g                    3D ray view: see the geodesics this camera casts\n"
"  h / p / v            hero shot, 4K / 720p / VGA -> kerr-hero-*.mp4\n"
"  r                    reset the camera\n"
"  s                    save a PNG snapshot\n"
"  q                    quit\n"
"\n"
"The disk does not turn in live mode: it changes faster than the tracer can\n"
"converge.  SPACE bakes a sequence instead -- the camera holds still, each\n"
"frame is taken to --bake-spp samples before the clock advances -- and loops\n"
"it back when you press SPACE again.  The window title always says what SPACE\n"
"will do next.  The finished sequence is written as an uncompressed .avi;\n"
"--video out.mp4 encodes H.265 instead, which is far smaller.\n"
"\n"
"The bake covers one ISCO orbit and closes exactly on itself; --no-loop keeps\n"
"the untouched pattern instead, which leaves a visible jump at the wrap.\n\n");
}

int main(int argc, char** argv) {
    Shared S;
    Scene& sc = S.scene;

    int  winW = 1280, winH = 720;
    Real diskIn = -1, diskOut = 18.0;
    Real targetFps = 30.0;
    int   bakeSpp = 64;      // samples per pixel before a baked frame is kept
    int   bakeMax = 180;     // frame cap; about 500 MB at 720p
    Real  bakeStep = -1;     // M between frames; <0 = derived from orbits/frames
    Real  bakeOrbits = 1.0;  // bake exactly this many ISCO orbits, then stop
    bool  bakeLoop   = true;  // --no-loop renders the true pattern with a visible seam
    int   bakeFrames = 96;   // frames per orbit
    std::string videoOut = "kerrbake.avi";   // .mp4 or .mkv here encodes H.265
    int    bakeCrf = 18;
    double heroSeconds = 60.0;   // H and V render this much footage
    int    heroSpp = 96;         // convergence per hero frame
    Real   heroSpan = 48.0;      // coordinate time the move covers, ~2 ISCO orbits
    double autoQuit = 0.0;      // scripted run: render for N seconds, snapshot, exit
    bool   autoOrbit = false;   // scripted camera motion, to exercise the moving path
    bool   autoRays  = false;   // scripted: open the ray view straight away
    std::string autoHero;       // scripted: "4k" or "vga", render and quit
    int    autoBake = 0;        // scripted: bake N frames, play, then quit
    double userExposure = 1.0, key = 1.3;
    int threads = 0;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto nextF = [&]() -> Real { return (i + 1 < argc) ? std::atof(argv[++i]) : 0; };
        auto nextI = [&]() -> int  { return (i + 1 < argc) ? std::atoi(argv[++i]) : 0; };
        if      (a == "--width")      winW = nextI();
        else if (a == "--height")     winH = nextI();
        else if (a == "--spin")       sc.kerr = Kerr(nextF());
        else if (a == "--dist")       S.cam.camR = nextF();
        else if (a == "--inc")        S.cam.incDeg = nextF();
        else if (a == "--fov")        S.cam.fovDeg = nextF();
        else if (a == "--rin")        diskIn = nextF();
        else if (a == "--rout")       diskOut = nextF();
        else if (a == "--tpeak")      sc.disk.tPeak = nextF();
        else if (a == "--albedo")     sc.disk.albedo = nextF();
        else if (a == "--edge")       sc.disk.edgeWidth = nextF();
        else if (a == "--turbulence") sc.disk.turbulence = nextF();
        else if (a == "--bounces")    sc.maxBounces = nextI();
        else if (a == "--sky-gain")   sc.sky.gain = nextF();
        else if (a == "--nostars")    sc.sky.enabled = false;
        else if (a == "--exposure")   userExposure = nextF();
        else if (a == "--fps")        targetFps = nextF();
        else if (a == "--time")       S.tObs = nextF();
        else if (a == "--bake-spp")   bakeSpp = nextI();
        else if (a == "--bake-step")  bakeStep = nextF();
        else if (a == "--bake-max")   bakeMax = nextI();
        else if (a == "--bake-orbits") bakeOrbits = nextF();
        else if (a == "--no-loop")    bakeLoop = false;
        else if (a == "--bake-frames") bakeFrames = nextI();
        else if (a == "--video")      videoOut = argv[++i];
        else if (a == "--no-video")   videoOut.clear();
        else if (a == "--crf")        bakeCrf = nextI();
        else if (a == "--hero-seconds") heroSeconds = nextF();
        else if (a == "--hero-spp")     heroSpp = nextI();
        else if (a == "--hero-span")    heroSpan = nextF();
        else if (a == "--autobake")   autoBake = nextI();
        else if (a == "--threads")    threads = nextI();
        else if (a == "--autoquit")   autoQuit = nextF();
        else if (a == "--autoorbit")  autoOrbit = true;
        else if (a == "--autorays")   autoRays = true;
        else if (a == "--autohero")   autoHero = argv[++i];
        else if (a == "--rtol")       sc.prop.rtol = nextF();
        else if (a == "--help" || a == "-h") {
            std::printf(
"kerrview -- interactive progressive viewer for the Kerr path tracer\n\n"
"  --width N --height N   window size            (1280 720)\n"
"  --spin A               spin a/M               (0.94)\n"
"  --dist R --inc DEG --fov DEG                  (40, 80, 40)\n"
"  --rin R --rout R --tpeak K --albedo A --edge W --turbulence F\n"
"  --bounces N --sky-gain F --nostars\n"
"  --exposure F --fps F --threads N --rtol F\n"
"  --autorays             open the 3D ray view at startup\n"
"  --autohero 4k|720p|vga render the hero shot and exit\n"
"  --hero-seconds F --hero-spp N --hero-span M\n"
"  --video F --no-video --crf N    bake output; .mp4/.mkv encode H.265\n");
            printControls();
            return 0;
        }
    }

    sc.disk.init(sc.kerr, diskIn, diskOut);
    sc.sky.seed = uint32_t(sc.seed * 40503u + 991u);
    sc.disk.noiseSeed = uint32_t(sc.seed * 2654435761u + 17u);
    sc.lamScale = spec::LAMBDA_SPAN / spec::ybarIntegral();

    sc.prop.kerr = &sc.kerr;
    sc.prop.disk = &sc.disk;
    sc.prop.rEscape = std::max<Real>(2000.0, S.cam.camR * 8);
    sc.prop.rCapture = std::max(sc.kerr.horizon() * 1.0005,
                                std::min(sc.kerr.photonOrbit(true), sc.disk.rIn));
    sc.prop.maxSteps = 60000;
    if (!(sc.prop.rtol > 0)) sc.prop.rtol = 1e-6;

    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        std::fprintf(stderr, "SDL_Init: %s\n", SDL_GetError());
        return 1;
    }
    SDL_Window* win = SDL_CreateWindow("kerrview",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, winW, winH,
        SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
    SDL_Renderer* ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);
    // Falling back to whatever renderer exists costs nothing here and is the
    // difference between running and not on a machine with no GPU path: X11
    // without GLX, a VM with 3D disabled, a plain SDL dummy driver.  There is
    // no GPU work to lose -- the frame is a CPU-side buffer that gets uploaded
    // once per present either way.  On Windows this never triggers, because
    // Direct3D is always there to be found.
    if (win && !ren) ren = SDL_CreateRenderer(win, -1, 0);
    if (!win || !ren) {
        std::fprintf(stderr, "SDL window/renderer: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }
    SDL_Texture* tex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_RGB24,
                                         SDL_TEXTUREACCESS_STREAMING, winW, winH);

    // One logical core is left to the UI thread: with every core saturated by
    // workers the event loop is starved and interaction feels worse than the
    // extra sample rate is worth.
    int hw = int(std::thread::hardware_concurrency());
    int nThreads = threads > 0 ? threads : std::max(1, hw - 1);
    if (nThreads < 1) nThreads = 1;
    S.nThreads = nThreads;

    S.rt.winW = winW; S.rt.winH = winH; S.rt.scale = 8;
    S.cam.aspect = Real(winW) / winH;
    S.display.assign(size_t(winW) * winH * 3, 0u);
    S.present = S.display;
    resetAccumulation(S);

    std::printf("kerrview -- %d x %d, %d worker threads\n", winW, winH, nThreads);
    std::printf("  spin a/M %.3f  horizon %.3f M  ISCO %.3f M  disk %.2f -> %.1f M\n",
                sc.kerr.a, sc.kerr.horizon(), sc.kerr.isco(true), sc.disk.rIn, sc.disk.rOut);
    printControls();
    std::fflush(stdout);

    std::vector<std::thread> pool;
    pool.reserve(nThreads);
    for (int i = 0; i < nThreads; ++i) pool.emplace_back(workerLoop, std::ref(S));

    // The main thread owns the authoritative camera and clock.  S.cam and
    // S.tObs are published from these only inside reconfigure(), i.e. with
    // every worker parked -- previously they were mutated in the event handler
    // while workers were reading them, which was a genuine data race.
    CameraParams camWanted = S.cam;
    Real         tWanted   = S.tObs;
    bool         sceneDirty = false;   // a change is pending publication

    Mode   mode = Mode::Live;
    RayVis rayVis;
    Hero   hero;
    // Leaving flies the inspection camera back to the viewpoint the rays were
    // launched from, and only then hands over to the traced image -- which is
    // rendered from exactly that pose, so the cut lands on a matching frame.
    bool              raysLeaving = false;
    Clock::time_point raysLeaveAt{};
    CameraParams      raysLeaveFrom{};
    const double      raysLeaveSecs = 1.1;
    Bake bake;
    bake.targetSpp = bakeSpp;
    bake.videoPath = videoOut;
    bake.crf       = bakeCrf;
    bake.orbitPeriod = TWO_PI / sc.disk.orbitOmega(sc.kerr.isco(true));
    {
        // Default to exactly one ISCO orbit, and tell the disk that is the loop
        // length.  Differential rotation means no span returns every ring to
        // its start on its own -- the outer disk needs twenty times as long as
        // the ISCO -- so the pattern is cross-faded against a copy one loop
        // older instead, which makes the last frame run back into the first
        // without touching any rotation rate.  See Disk::mottle.
        int wanted = std::max(2, int(bakeFrames * bakeOrbits + 0.5));
        bake.step = (bakeStep > 0) ? bakeStep : (bake.orbitPeriod * bakeOrbits) / wanted;
        bake.maxFrames = std::min(bakeMax, wanted);
        bake.loopSpan = bake.step * bake.maxFrames;
    }
    const CameraParams home = camWanted;

    bool  dragOrbit = false, dragPan = false;
    auto  lastInput = Clock::now() - std::chrono::seconds(10);
    auto  lastStats = Clock::now();
    auto  lastTitle = Clock::now();
    auto  lastStatus = Clock::now();
    uint64_t statSamples = 0;
    double raysPerSec = 2.0e6;          // seeded, corrected after the first second
    // Meter once, coarsely and synchronously, before anything is shown.
    double logExposure = meterInitial(sc, S.cam, key);
    S.exposure.store(std::exp(logExposure) * userExposure, std::memory_order_relaxed);
    double fps = 0.0;
    auto  lastFrame = Clock::now();
    auto  lastPresent = Clock::now() - std::chrono::seconds(1);
    auto  lastPlayTick = Clock::now();
    bool  running = true;
    int   snapshotIndex = 0;

    // A scripted hero shot has nothing to do afterwards, so it exits.
    // A scripted hero exits when it is done, unless --autoquit says the run
    // has something else to do afterwards -- which is how the return to
    // ordinary viewing gets exercised.
    const bool heroThenQuit = !autoHero.empty() && autoQuit <= 0;

    // Put the viewer back the way it was, and close the container.
    auto heroFinish = [&](const char* why) {
        if (hero.videoOpen) {
            if (!hero.vw.close())
                std::printf("\n  hero: the encoder reported a failure; %s may be truncated\n",
                            hero.path.c_str());
            hero.videoOpen = false;
        }
        reconfigure(S, [&] {
            S.rt.winW = hero.saveW; S.rt.winH = hero.saveH; S.rt.scale = hero.saveScale;
            S.display.assign(size_t(hero.saveW) * hero.saveH * 3, 0u);
            S.present = S.display;
            camWanted = hero.saveCam;
            tWanted   = hero.saveT;
            S.cam        = camWanted;
            S.cam.aspect = Real(hero.saveW) / std::max(1, hero.saveH);
            S.tObs       = tWanted;
            // passLimit is normally maintained by the regime block further
            // down, which this mode never reaches.  Leaving a stale cap here
            // parks every worker the moment the counter passes it.
            S.passLimit.store(0, std::memory_order_relaxed);
            resetAccumulation(S);
        });
        std::printf("\n  hero: %s -- %d frames to %s (%.0f s)\n",
                    why, hero.idx, hero.path.c_str(), since(hero.began));
        std::fflush(stdout);
        hero.active = false;
        mode = Mode::Live;
        // Deliberately not touching lastInput.  Doing so makes `moving` true,
        // the regime block then caps dispatch to one item per tile, and once
        // moving lapses nothing fires again to lift the cap -- the workers
        // finish a single pass and park for good.  Marking the scene dirty
        // instead gets one settled pass through that block, which clears it.
        sceneDirty = true;
        if (heroThenQuit) running = false;
    };

    // Set the move going at the given resolution.  The estimate is printed
    // before the first frame on purpose: a minute of 4K is a very long job, and
    // ESC is a lot cheaper to press now than in an hour.
    auto heroStart = [&](int hw, int hh, const char* label) {
        if (mode != Mode::Live) return;
        hero = Hero{};
        hero.w = hw; hero.h = hh;
        hero.fps = 24.0;
        hero.frames = std::max(2, int(heroSeconds * hero.fps + 0.5));
        hero.spp = std::max(1, heroSpp);
        hero.tSpan = heroSpan;
        hero.t0 = tWanted;
        hero.path = std::string("kerr-hero-") + label + ".mp4";
        hero.saveW = S.rt.winW; hero.saveH = S.rt.winH; hero.saveScale = S.rt.scale;
        hero.saveCam = camWanted; hero.saveT = tWanted;
        hero.began = Clock::now();

        hero.vw.crf = bakeCrf;
        hero.videoOpen = hero.vw.open(hero.path, hw, hh, hero.fps);
        if (!hero.videoOpen) {
            std::printf("  hero: cannot open %s for writing\n", hero.path.c_str());
            std::fflush(stdout);
            return;
        }

        double rays = double(hw) * hh * hero.spp * hero.frames;
        double est  = rays / std::max(1.0e5, raysPerSec);
        std::printf("\n  hero: %dx%d, %d frames at %.0f fps (%.0f s of footage), %d spp\n"
                    "        %.1f Gray total, about %.1f h at the current rate -- ESC to abandon\n",
                    hw, hh, hero.frames, hero.fps, hero.frames / hero.fps, hero.spp,
                    rays / 1e9, est / 3600.0);
        std::fflush(stdout);

        camWanted = heroPose(home, 0.0);
        tWanted   = hero.t0;
        reconfigure(S, [&] {
            S.rt.winW = hw; S.rt.winH = hh; S.rt.scale = 1;
            S.display.assign(size_t(hw) * hh * 3, 0u);
            S.present = S.display;
            S.cam        = camWanted;
            S.cam.aspect = Real(hw) / hh;
            S.tObs       = tWanted;
            S.scene.disk.loopPeriod = 0;   // no loop trick: the disk simply turns
            S.passLimit.store(0, std::memory_order_relaxed);   // converge, uncapped
            resetAccumulation(S);
        });
        hero.active = true;
        mode = Mode::Hero;
    };

    auto startTime = Clock::now();
    while (running) {
        // Meter first, on whatever the workers produced since the last reset.
        // Doing this after the reset below always reads an empty buffer, which
        // is exactly how the exposure used to end up meaningless while moving.
        logExposure = updateExposure(S, key, userExposure, logExposure);

        // Regime for this frame, decided before events so the snapshot paths
        // below can tell which buffer is on screen.
        bool moving = (mode == Mode::Live) && since(lastInput) < 0.25;
        bool camChanged = false;

        // Advancing the disk invalidates the accumulated image exactly as a
        // camera move does, so it goes through the same reset path.  That also
        // means an animating disk never settles -- which is not a limitation
        // but arithmetic: a changing scene has nothing to accumulate.  Pause
        // with space to let it converge.
        // The viewer deliberately does not animate the disk.
        //
        // A rotating disk changes faster than the path tracer can converge, so
        // an animating viewer is permanently stuck at one sample per pixel and
        // can never show a settled image -- progressive refinement and
        // animation are directly at odds.  The scene here is therefore a single
        // instant (`--time`), which lets the image converge whenever you stop
        // moving the camera.  Rotation lives in the batch renderer, where every
        // frame of a sequence can be taken to full convergence: see
        // `kerr.exe --frames`.

        if (mode == Mode::Rays && raysLeaving) {
            double u = std::min(1.0, since(raysLeaveAt) / raysLeaveSecs);
            double e = u * u * (3 - 2 * u);                  // smoothstep
            camWanted = lerpPose(raysLeaveFrom, rayVis.source, Real(e));
            if (u >= 1.0) {
                camWanted   = rayVis.source;
                mode        = Mode::Live;
                raysLeaving = false;
                camChanged  = true;
                lastInput   = Clock::now();
                std::printf("  rays: back to live\n");
                std::fflush(stdout);
            }
        }

        // --- the hero shot ------------------------------------------------
        if (mode == Mode::Hero && hero.active) {
            // The exposure is deliberately *not* frozen here, unlike a bake.
            // A bake holds the camera still, so freezing keeps playback from
            // flickering; this move triples the disk's share of the frame on the
            // way in, and a level metered from the wide shot burns the close-up
            // to white.  The viewer's metering is already smoothed over time,
            // which is what keeps the adjustment from showing.
            if (S.minCount.load(std::memory_order_relaxed) >= uint32_t(hero.spp)) {
                bool finished = false;
                // The frame is taken with the workers parked; reading display
                // while they paint it would tear.
                reconfigure(S, [&] {
                    if (hero.videoOpen) hero.vw.addFrame(S.display.data());
                    ++hero.idx;
                    if (hero.idx >= hero.frames) { finished = true; return; }
                    double u   = double(hero.idx) / double(std::max(1, hero.frames - 1));
                    camWanted  = heroPose(home, u);
                    tWanted    = hero.t0 + hero.tSpan * (double(hero.idx) / hero.frames);
                    S.cam        = camWanted;
                    S.cam.aspect = Real(hero.w) / hero.h;
                    S.tObs       = tWanted;
                    S.passLimit.store(0, std::memory_order_relaxed);
                    resetAccumulation(S);
                });
                if (hero.idx % 10 == 0 || finished) {
                    double el = since(hero.began);
                    double eta = hero.idx > 0 ? el / hero.idx * (hero.frames - hero.idx) : 0;
                    std::printf("\r  hero  %d/%d   %.0f%%   %.0f s elapsed, %.0f s left    ",
                                hero.idx, hero.frames, 100.0 * hero.idx / hero.frames, el, eta);
                    std::fflush(stdout);
                }
                if (finished) heroFinish("done");
            }
        }

        // Scripted hero shot, so it can be rendered without a hand on the keyboard.
        if (!autoHero.empty() && mode == Mode::Live && !hero.active
            && since(startTime) > 0.4) {
            if      (autoHero == "4k")   heroStart(3840, 2160, "4k");
            else if (autoHero == "720p") heroStart(1280, 720,  "720p");
            else                         heroStart(640,  480,  "vga");
            autoHero.clear();
        }

        // Scripted ray view, so it can be exercised without a hand on the mouse.
        if (autoRays && mode == Mode::Live && !rayVis.captured && since(startTime) > 0.4) {
            reconfigure(S, [&] {
                captureRays(rayVis, sc, camWanted);
                rayVis.imgPixels = S.display;
                rayVis.imgW = winW;
                rayVis.imgH = winH;
            });
            rebuildRayTexture(ren, rayVis);
            mode = Mode::Rays;
            std::printf("  [auto] ray view: %zu trajectories\n", rayVis.lines.size());
            std::fflush(stdout);
        }

        // Scripted bake, so the workflow can be exercised headlessly.
        if (autoBake > 0) {
            if (mode == Mode::Live && since(startTime) > 0.4) {
                bake.frames.clear();
                bake.t0 = tWanted;
                bake.lockedExposure = S.exposure.load();
                mode = Mode::Baking;
                sceneDirty = true;
                std::printf("  [auto] baking %d frames\n", autoBake);
                std::fflush(stdout);
            } else if (mode == Mode::Baking && int(bake.frames.size()) >= autoBake) {
                mode = Mode::Playing;
                bake.playHead = 0;
                for (size_t i = 0; i < bake.frames.size(); ++i) {
                    char nm[128];
                    std::snprintf(nm, sizeof nm, "kerrbake-%04zu.png", i);
                    img::writePng(nm, winW, winH, bake.frames[i]);
                }
                writeBakeVideo(bake, winW, winH);
                std::printf("  [auto] playing %zu frames (exported for inspection)\n",
                            bake.frames.size());
                std::fflush(stdout);
            }
        }

        // Scripted motion, so the interactive path can be exercised headlessly.
        if (autoOrbit && since(startTime) < autoQuit * 0.6) {
            camWanted.phiDeg += 1.5;
            camChanged = true;
        }
        if (autoQuit > 0 && since(startTime) > autoQuit) {
            std::vector<uint8_t> copy;
            if (mode == Mode::Playing && !bake.frames.empty())
                copy = bake.frames[std::min(bake.frames.size() - 1, size_t(bake.playHead))];
            else
                reconfigure(S, [&] { copy = moving ? S.present : S.display; });
            if (mode == Mode::Rays) {
                savePresented(ren, winW, winH, "kerrview-rays.png");
                std::printf("  autoquit: wrote kerrview-rays.png, %zu trajectories\n",
                            rayVis.lines.size());
            } else {
                img::writePng("kerrview-auto.png", winW, winH, copy);
                std::printf("  autoquit: wrote kerrview-auto.png at %u spp, scale 1/%d, %.2f Mrays/s\n",
                            S.minCount.load(), S.rt.scale, raysPerSec / 1e6);
            }
            running = false;
        }

        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            switch (e.type) {
            case SDL_QUIT: running = false; break;

            case SDL_WINDOWEVENT:
                if (e.window.event == SDL_WINDOWEVENT_SIZE_CHANGED) {
                    winW = e.window.data1; winH = e.window.data2;
                    SDL_DestroyTexture(tex);
                    tex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_RGB24,
                                            SDL_TEXTUREACCESS_STREAMING, winW, winH);
                    reconfigure(S, [&] {
                        S.rt.winW = winW; S.rt.winH = winH;
                        S.cam.aspect = Real(winW) / std::max(1, winH);
                        S.display.assign(size_t(winW) * winH * 3, 0u);
                        S.present = S.display;
                        camWanted.aspect = S.cam.aspect;
                        resetAccumulation(S);
                    });
                    lastInput = Clock::now();
                }
                break;

            case SDL_MOUSEBUTTONDOWN:
                if (e.button.button == SDL_BUTTON_LEFT)  dragOrbit = true;
                if (e.button.button == SDL_BUTTON_RIGHT) dragPan   = true;
                break;
            case SDL_MOUSEBUTTONUP:
                if (e.button.button == SDL_BUTTON_LEFT)  dragOrbit = false;
                if (e.button.button == SDL_BUTTON_RIGHT) dragPan   = false;
                break;

            case SDL_MOUSEMOTION:
                // Rays mode orbits the inspection camera; baking and playback
                // assume a fixed one.
                if (mode != Mode::Live && mode != Mode::Rays) break;
                if (mode == Mode::Rays && raysLeaving) break;
                if (dragOrbit && (e.motion.xrel || e.motion.yrel)) {
                    // Negated to match the other two axes: screen-right is the
                    // direction of increasing phi, so orbiting with +xrel swings
                    // the camera the same way the cursor went and the scene
                    // parallaxes against it.  Drag should carry the hole with it.
                    camWanted.phiDeg -= e.motion.xrel * 0.25;
                    camWanted.incDeg = clampf(camWanted.incDeg - e.motion.yrel * 0.25, 1.0, 179.0);
                    camChanged = true;
                } else if (dragPan && (e.motion.xrel || e.motion.yrel)) {
                    camWanted.yawDeg   = clampf(camWanted.yawDeg   - e.motion.xrel * 0.05, -80.0, 80.0);
                    camWanted.pitchDeg = clampf(camWanted.pitchDeg + e.motion.yrel * 0.05, -80.0, 80.0);
                    camChanged = true;
                }
                break;

            case SDL_MOUSEWHEEL: {
                if (mode != Mode::Live && mode != Mode::Rays) break;
                bool ctrl = (SDL_GetModState() & KMOD_CTRL) != 0;
                if (ctrl) {
                    camWanted.fovDeg = clampf(camWanted.fovDeg * std::exp(-e.wheel.y * 0.08), 2.0, 120.0);
                } else {
                    Real minR = sc.kerr.horizon() * 2.0 + 1.0;
                    camWanted.camR = clampf(camWanted.camR * std::exp(-e.wheel.y * 0.10), minR, 4000.0);
                    sc.prop.rEscape = std::max<Real>(2000.0, camWanted.camR * 8);
                }
                camChanged = true;
                break;
            }

            case SDL_KEYDOWN:
                switch (e.key.keysym.sym) {
                case SDLK_q: running = false; break;

                // One key for the whole workflow.  Whatever is on screen, space
                // does the next sensible thing, and the title says what that is.
                case SDLK_SPACE:
                    if (mode == Mode::Live) {
                        bake.frames.clear();
                        bake.t0 = tWanted;
                        bake.lockedExposure = S.exposure.load();
                        mode = Mode::Baking;
                        sceneDirty = true;
                        std::printf("  BAKING: %d spp per frame, %.2f M apart, up to %d frames."
                                    "  SPACE to stop and play.\n",
                                    bake.targetSpp, double(bake.step), bake.maxFrames);
                    } else if (mode == Mode::Baking) {
                        if (bake.frames.size() >= 2) {
                            mode = Mode::Playing;
                            bake.playHead = 0;
                            writeBakeVideo(bake, winW, winH);
                            std::printf("  PLAYING %zu frames at %.0f fps (%.1f M, %.2f ISCO orbits)."
                                        "  SPACE for live.\n",
                                        bake.frames.size(), bake.fps, double(bake.span()),
                                        double(bake.span() / (TWO_PI / sc.disk.orbitOmega(sc.kerr.isco(true)))));
                        } else {
                            std::printf("  need at least 2 frames to play; still baking\n");
                        }
                    } else {
                        mode = Mode::Live;
                        tWanted = bake.t0;
                        sceneDirty = true;
                        std::printf("  live\n");
                    }
                    std::fflush(stdout);
                    lastInput = Clock::now();
                    break;

                // Escape backs out of whatever is happening, without quitting.
                case SDLK_ESCAPE:
                    if (mode == Mode::Hero && hero.active) { heroFinish("abandoned"); break; }
                    if (mode == Mode::Rays) {
                        mode = Mode::Live;
                        raysLeaving = false;
                        camWanted = rayVis.source;
                        camChanged = true;
                        std::printf("  rays: back to live\n");
                        std::fflush(stdout);
                        lastInput = Clock::now();
                        break;
                    }
                    if (mode != Mode::Live) {
                        mode = Mode::Live;
                        tWanted = bake.t0;
                        sceneDirty = true;
                        std::printf("  cancelled, back to live\n");
                        std::fflush(stdout);
                        lastInput = Clock::now();
                    }
                    break;
                case SDLK_r: camWanted = home; camChanged = true; break;
                case SDLK_h: heroStart(3840, 2160, "4k");   break;
                case SDLK_p: heroStart(1280, 720,  "720p"); break;
                case SDLK_v: heroStart(640,  480,  "vga");  break;
                case SDLK_g:
                    if (mode == Mode::Live) {
                        // Capture with the workers parked: the scene is only
                        // safe to read while nothing is tracing it.
                        reconfigure(S, [&] {
                            captureRays(rayVis, sc, camWanted);
                            rayVis.imgPixels = S.display;    // the frame these rays made
                            rayVis.imgW = winW;
                            rayVis.imgH = winH;
                        });
                        rebuildRayTexture(ren, rayVis);
                        mode = Mode::Rays;
                        std::printf("  rays: %zu trajectories from r = %.1f M, inc %.0f deg"
                                    "  (drag to orbit, [ / ] density, G or ESC to return)\n",
                                    rayVis.lines.size(), double(camWanted.camR),
                                    double(camWanted.incDeg));
                        std::fflush(stdout);
                    } else if (mode == Mode::Rays && !raysLeaving) {
                        raysLeaving   = true;
                        raysLeaveAt   = Clock::now();
                        raysLeaveFrom = camWanted;
                        std::printf("  rays: flying back\n");
                        std::fflush(stdout);
                    }
                    break;
                case SDLK_LEFTBRACKET:
                case SDLK_RIGHTBRACKET: {
                    if (mode == Mode::Rays) {
                        int d = (e.key.keysym.sym == SDLK_RIGHTBRACKET) ? 2 : -2;
                        rayVis.gridX = std::max(3, std::min(41, rayVis.gridX + d));
                        rayVis.gridY = std::max(3, std::min(23, rayVis.gridY + (d > 0 ? 1 : -1)));
                        reconfigure(S, [&] { captureRays(rayVis, sc, rayVis.source); });
                        std::printf("  rays: %d x %d = %zu trajectories\n",
                                    rayVis.gridX, rayVis.gridY, rayVis.lines.size());
                        std::fflush(stdout);
                        break;
                    }
                    if (mode != Mode::Live) break;
                    Real d = (e.key.keysym.sym == SDLK_RIGHTBRACKET) ? 0.02 : -0.02;
                    reconfigure(S, [&] {
                        sc.kerr = Kerr(clampf(sc.kerr.a + d, -0.999, 0.999));
                        sc.disk.init(sc.kerr, diskIn, diskOut);
                        sc.prop.rCapture = std::max(sc.kerr.horizon() * 1.0005,
                                                    std::min(sc.kerr.photonOrbit(true), sc.disk.rIn));
                        resetAccumulation(S);
                    });
                    std::printf("  spin a/M = %.3f   ISCO %.3f M\n", sc.kerr.a, sc.kerr.isco(true));
                    std::fflush(stdout);
                    lastInput = Clock::now();
                    break;
                }
                case SDLK_COMMA:
                case SDLK_PERIOD: {
                    Real f = (e.key.keysym.sym == SDLK_PERIOD) ? 1.1 : 1.0 / 1.1;
                    if (mode == Mode::Playing) {   // same keys, obvious meaning
                        bake.fps = clampf(bake.fps * f, 1.0, 120.0);
                        std::printf("  playback %.0f fps\n", bake.fps);
                        std::fflush(stdout);
                        break;
                    }
                    if (mode != Mode::Live) break;
                    diskOut = clampf(diskOut * f, sc.disk.rIn * 1.2, 400.0);
                    reconfigure(S, [&] {
                        sc.disk.init(sc.kerr, diskIn, diskOut);
                        resetAccumulation(S);
                    });
                    std::printf("  disk outer radius = %.1f M\n", diskOut);
                    std::fflush(stdout);
                    lastInput = Clock::now();
                    break;
                }
                case SDLK_MINUS:  userExposure *= 1.0 / 1.25; break;
                case SDLK_EQUALS: userExposure *= 1.25;       break;
                case SDLK_b: {
                    if (mode != Mode::Live) break;
                    reconfigure(S, [&] {
                        sc.maxBounces = (sc.maxBounces + 1) % 5;
                        resetAccumulation(S);
                    });
                    std::printf("  scattering bounces = %d\n", sc.maxBounces);
                    std::fflush(stdout);
                    lastInput = Clock::now();
                    break;
                }
                case SDLK_s: {
                    if (mode == Mode::Rays) {
                        if (savePresented(ren, winW, winH, "kerrview-rays.png"))
                            std::printf("  wrote kerrview-rays.png\n");
                        else
                            std::printf("  could not read the framebuffer back: %s\n", SDL_GetError());
                        std::fflush(stdout);
                        break;
                    }
                    // During playback, save the whole baked sequence -- having
                    // just watched it, that is obviously what "save" means.
                    if (mode == Mode::Playing && !bake.frames.empty()) {
                        int written = 0;
                        for (size_t i = 0; i < bake.frames.size(); ++i) {
                            char name[128];
                            std::snprintf(name, sizeof name, "kerrbake-%04zu.png", i);
                            if (img::writePng(name, winW, winH, bake.frames[i])) ++written;
                        }
                        std::printf("  wrote %d frames as kerrbake-0000.png ...\n", written);
                        std::fflush(stdout);
                        break;
                    }
                    char name[128];
                    std::snprintf(name, sizeof name, "kerrview-%03d.png", snapshotIndex++);
                    std::vector<uint8_t> copy;
                    reconfigure(S, [&] { copy = moving ? S.present : S.display; });
                    if (img::writePng(name, winW, winH, copy))
                        std::printf("  wrote %s  (%u spp)\n", name, S.minCount.load());
                    else
                        std::printf("  failed to write %s\n", name);
                    std::fflush(stdout);
                    break;
                }
                default: break;
                }
                break;
            }
        }
        if (camChanged) lastInput = Clock::now();

        // --- bake and playback ---------------------------------------------
        //
        // Baking forces the settled regime: full resolution, unlimited
        // accumulation, camera held still.  When the frame reaches its sample
        // target it is captured and the clock steps on.
        if (mode == Mode::Baking) {
            // Freeze the exposure for the whole sequence.  Re-metering per
            // frame would let the level drift between them, which reads as
            // flicker during playback.
            S.exposure.store(bake.lockedExposure, std::memory_order_relaxed);

            if (S.minCount.load() >= uint32_t(bake.targetSpp)) {
                bake.frames.push_back(S.display);
                tWanted += bake.step;
                sceneDirty = true;

                Real period = TWO_PI / sc.disk.orbitOmega(sc.kerr.isco(true));
                Real orbits = bake.span() / period;
                std::printf("  frame %3zu   t = %7.2f M   %.3f ISCO orbits   %.0f MB%s\n",
                            bake.frames.size(), double(tWanted - bake.t0), double(orbits),
                            bake.bytes() / 1.0e6,
                            (std::fabs(orbits - std::round(orbits)) < 0.5 * double(bake.step) / double(period)
                             && orbits >= 0.9) ? "   <- whole orbit, good loop point" : "");
                std::fflush(stdout);

                if (int(bake.frames.size()) >= bake.maxFrames) {
                    mode = Mode::Playing;
                    bake.playHead = 0;
                    writeBakeVideo(bake, winW, winH);
                    std::printf("  frame limit reached; PLAYING %zu frames\n", bake.frames.size());
                    std::fflush(stdout);
                }
            }
        }

        // Playback needs no rays at all, so park the workers rather than leave
        // 31 threads spinning on an image nobody is looking at.
        {
            bool wantIdle = (mode == Mode::Playing || mode == Mode::Rays);
            std::unique_lock<std::mutex> lk(S.mtx);
            if (S.idle != wantIdle) { S.idle = wantIdle; lk.unlock(); S.cv.notify_all(); }
        }

        if (mode == Mode::Playing && !bake.frames.empty()) {
            bake.playHead += bake.fps * std::min(0.25, since(lastPlayTick));
            lastPlayTick = Clock::now();
            if (bake.playHead >= double(bake.frames.size())) bake.playHead = 0;
        } else {
            lastPlayTick = Clock::now();
        }

        // --- choose the regime -------------------------------------------
        moving = (mode == Mode::Live) && since(lastInput) < 0.25;
        int desiredScale;
        if (moving) {
            // Pick the coarsest-to-finest scale that still fits one sample per
            // pixel into the frame budget at the measured throughput.
            double affordable = std::max(1000.0, raysPerSec / std::max(1.0, targetFps));
            double ratio = double(winW) * winH / affordable;
            desiredScale = std::max(1, int(std::ceil(std::sqrt(std::max(1.0, ratio)))));
            desiredScale = std::min(desiredScale, 16);
        } else {
            desiredScale = 1;
        }

        // Coalesce scene changes until the current one has actually been drawn.
        //
        // Resetting on every UI frame is what made an animating disk look
        // stalled: the display loop runs at 60 Hz while a full pass completes
        // at ~40 Hz, so the accumulator was wiped before any pass finished and
        // the frame was permanently about half current and half leftovers.
        // Waiting for a completed pass means every image shown is a complete
        // render of one scene state, at the cost of at most one pass of input
        // latency (~25 ms).  A scale change still applies at once, since the
        // buffers have to be resized anyway.
        bool passComplete = S.work.load(std::memory_order_relaxed) >= uint64_t(S.numTiles);
        if (camChanged && mode != Mode::Rays && mode != Mode::Hero) sceneDirty = true;

        if ((passComplete && (sceneDirty || moving)) || desiredScale != S.rt.scale) {
            reconfigure(S, [&] {
                // Publish the finished pass.  While moving, the texture is
                // uploaded from this copy rather than from the live buffer, so
                // every frame shown has exactly one sample per pixel at one
                // instant.  Reading the live buffer instead shows a mix of two
                // passes, which appears as rectangular tile patches of
                // differing noise and differing epoch.
                S.present = S.display;
                if (sceneDirty) {
                    S.rt.scale = desiredScale;
                    S.cam  = camWanted;   // published with every worker parked
                    S.tObs = tWanted;
                    // The loop cross-fade belongs to a sequence, so it is only
                    // in force while one is being captured.  Live mode shows
                    // the plain advected pattern at whatever instant it is on.
                    S.scene.disk.loopPeriod =
                        (bakeLoop && mode == Mode::Baking) ? bake.loopSpan : Real(0);
                    resetAccumulation(S);
                }
                // Exactly one sample per tile while moving; unlimited once
                // settled, so the image goes on converging.
                S.passLimit.store(moving ? uint64_t(S.numTiles) : 0,
                                  std::memory_order_relaxed);
            });
            sceneDirty = false;
        }
        if (desiredScale != S.rt.scale) {
            reconfigure(S, [&] {
                S.rt.scale = desiredScale;
                resetAccumulation(S);
                S.passLimit.store(moving ? uint64_t(S.numTiles) : 0, std::memory_order_relaxed);
            });
        }

        // --- throughput and exposure --------------------------------------
        double dtStats = since(lastStats);
        if (dtStats > 0.25) {
            uint64_t now = S.samplesTraced.load(std::memory_order_relaxed);
            double inst = double(now - statSamples) / dtStats;
            raysPerSec += (inst - raysPerSec) * 0.35;
            statSamples = now;
            lastStats = Clock::now();
        }
        // Cheapest useful progress figure: samples completed per render pixel.
        {
            uint64_t w = S.work.load(std::memory_order_relaxed);
            S.minCount.store(uint32_t(w / uint64_t(std::max(1, S.numTiles))), std::memory_order_relaxed);
        }

        // --- present -------------------------------------------------------
        //
        // The loop deliberately does not block on vsync.  A capped pass can
        // finish part way through a display interval, and if publication waits
        // for the next vsync every worker sits idle until then -- measured at a
        // third of total throughput, which the adaptive scale then pays for in
        // resolution.  Polling fast and presenting on a timer instead keeps
        // both: the pass is published as soon as it completes, while the screen
        // still updates at 60 Hz.
        if (since(lastPresent) >= 1.0 / 60.0) {
            // Playback shows a baked frame; otherwise the last complete pass
            // while moving, or the live buffer while settled so refinement is
            // visible as it happens.
            // The visualiser draws vectors rather than a frame buffer, so it
            // bypasses the streaming texture entirely.
            if (mode == Mode::Hero) {
                // The render target is the hero resolution, not the window's, so
                // there is nothing here the streaming texture could show.  A bar
                // is more use than a stretched corner of a 4K frame.
                SDL_SetRenderDrawColor(ren, 8, 8, 10, 255);
                SDL_RenderClear(ren);
                double f = double(hero.idx) / std::max(1, hero.frames);
                int bw = winW - 80, bx = 40, by = winH / 2 - 14;
                SDL_Rect bg{bx, by, bw, 28};
                SDL_SetRenderDrawColor(ren, 38, 38, 46, 255);
                SDL_RenderFillRect(ren, &bg);
                SDL_Rect fg{bx, by, int(bw * f), 28};
                SDL_SetRenderDrawColor(ren, 255, 176, 64, 255);
                SDL_RenderFillRect(ren, &fg);
            } else if (mode == Mode::Rays) {
                drawRayVis(ren, rayVis, camWanted, winW, winH);
            } else {
                const uint8_t* src;
                if (mode == Mode::Playing && !bake.frames.empty()) {
                    size_t i = std::min(bake.frames.size() - 1, size_t(bake.playHead));
                    src = bake.frames[i].data();
                } else {
                    src = moving ? S.present.data() : S.display.data();
                }
                SDL_UpdateTexture(tex, nullptr, src, winW * 3);
                SDL_RenderCopy(ren, tex, nullptr, nullptr);
            }
            SDL_RenderPresent(ren);
            lastPresent = Clock::now();

            double dtFrame = since(lastFrame);
            lastFrame = Clock::now();
            if (dtFrame > 0) fps += (1.0 / dtFrame - fps) * 0.1;
        } else {
            // Nothing to draw yet; yield rather than spin a core on it.
            std::this_thread::sleep_for(std::chrono::microseconds(300));
        }

        if (autoQuit > 0 && since(lastStatus) > 0.5) {
            double passSamples = double(S.rt.rw) * S.rt.rh;
            double passHz = raysPerSec / std::max(1.0, passSamples);
            std::printf("  t=%4.1fs  %-8s scale 1/%-2d  render %dx%d  %2u spp  ui %.0f fps  "
                        "first-pass %.1f Hz  %.2f Mrays/s\n",
                        since(startTime), moving ? "moving" : "settled", S.rt.scale,
                        S.rt.rw, S.rt.rh, S.minCount.load(), fps, passHz, raysPerSec / 1e6);
            std::printf("           t = %7.2f M   exposure %.4g   metered from %.0f%% of the frame\n",
                        double(S.tObs), S.exposure.load(), 100.0 * S.meterCovered);
            std::fflush(stdout);
            lastStatus = Clock::now();
        }
        if (since(lastTitle) > 0.2) {
            // The title is the only text surface available, so it carries the
            // whole interface: what state we are in, how it is progressing, and
            // -- last, always -- what SPACE will do next.
            char title[320];
            Real period = TWO_PI / sc.disk.orbitOmega(sc.kerr.isco(true));
            switch (mode) {
            case Mode::Baking:
                std::snprintf(title, sizeof title,
                    "BAKING  |  frame %zu/%d  |  %u/%d spp  |  %.2f ISCO orbits captured  |  "
                    "%.0f MB  |  SPACE: stop and play    ESC: cancel",
                    bake.frames.size() + 1, bake.maxFrames,
                    S.minCount.load(), bake.targetSpp,
                    double(bake.span() / period), bake.bytes() / 1.0e6);
                break;
            case Mode::Hero:
                std::snprintf(title, sizeof title,
                    "HERO  |  %d/%d  |  %dx%d  |  %u/%d spp  |  %.0f%%  |  ESC to abandon",
                    hero.idx, hero.frames, hero.w, hero.h,
                    S.minCount.load(), hero.spp,
                    100.0 * hero.idx / std::max(1, hero.frames));
                break;
            case Mode::Rays:
                std::snprintf(title, sizeof title,
                    "RAYS  |  %zu trajectories  |  %d x %d  |  from r=%.0fM inc=%.0f  |  "
                    "drag to orbit  |  G: back to live",
                    rayVis.lines.size(), rayVis.gridX, rayVis.gridY,
                    double(rayVis.source.camR), double(rayVis.source.incDeg));
                break;
            case Mode::Playing:
                std::snprintf(title, sizeof title,
                    "PLAYING  |  frame %zu/%zu  |  %.0f fps  |  %.2f ISCO orbits  |  "
                    ", / . speed  |  SPACE: back to live",
                    size_t(bake.playHead) + 1, bake.frames.size(), bake.fps,
                    double(bake.span() / period));
                break;
            default:
                std::snprintf(title, sizeof title,
                    "kerrview  |  %s 1/%d res  |  %u spp  |  %.1f Mrays/s  |  "
                    "a=%.2f r=%.0fM inc=%.0f  |  SPACE: bake an animation",
                    moving ? "moving" : "settled", S.rt.scale, S.minCount.load(),
                    raysPerSec / 1e6,
                    sc.kerr.a, double(camWanted.camR), double(camWanted.incDeg));
                break;
            }
            SDL_SetWindowTitle(win, title);
            lastTitle = Clock::now();
        }
    }

    {
        std::lock_guard<std::mutex> lk(S.mtx);
        S.quit = true;
        S.paused = false;
    }
    S.cv.notify_all();
    for (auto& t : pool) t.join();

    if (rayVis.imgTex) SDL_DestroyTexture(rayVis.imgTex);
    SDL_DestroyTexture(tex);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}
