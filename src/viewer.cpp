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
    bool quit   = false;
    int  activeWorkers = 0;

    // Mutated only while paused and activeWorkers == 0.
    Scene        scene;
    CameraParams cam;
    Target       rt;

    std::vector<double>   accum;    // rw*rh*3, CIE XYZ
    std::vector<uint32_t> count;    // rw*rh
    std::vector<uint8_t>  display;  // winW*winH*3, persists across resets

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
    std::atomic<uint64_t> samplesTraced{0};
    std::atomic<double>   exposure{1.0};
    std::atomic<uint32_t> minCount{0};      // samples per pixel, for the title
    double meterCovered = 0, meterAnchor = 0, meterTarget = 0;
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
    auto shade = [&](size_t idx, int x, int y, Real lambda, Real L) {
        double* a = &S.accum[idx * 3];
        if (L > 0) {
            Vec3 c = spec::cieXYZ(lambda);
            a[0] += c.x * L; a[1] += c.y * L; a[2] += c.z * L;
        }
        uint32_t n = ++S.count[idx];

        // Tone map this pixel now, while its data is hot in cache.
        Real inv = sc.lamScale / n * expScale;
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
            Real     lam[LANES];
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
                lam[j] = spec::LAMBDA_MIN + spec::LAMBDA_SPAN * ul;

                Real sx = 2.0 * (px + u1) / rt.rw - 1.0;
                Real sy = 1.0 - 2.0 * (y + u2) / rt.rh;
                gp[j] = cam.ray(sx, sy);
            }

            Real rad[LANES];
            tracePacket(sc.kerr, sc.disk, sc.sky, sc.prop, gp, lam, rgs,
                        sc.maxBounces, n, rad);
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
            S.cv.wait(lk, [&] { return !S.paused || S.quit; });
            if (S.quit) return;
            ++S.activeWorkers;
        }
        // Safe to touch scene/buffers unlocked: the main thread only mutates
        // them once activeWorkers has fallen to zero, which cannot happen
        // while we are counted.
        uint64_t w = S.work.fetch_add(1, std::memory_order_relaxed);
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
        lum.push_back(S.accum[i * 3 + 1] / S.count[i]);
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
                Real lambda = spec::LAMBDA_MIN + spec::LAMBDA_SPAN * ul;
                Real sx = 2.0 * (x + u1) / w - 1.0;
                Real sy = 1.0 - 2.0 * (y + u2) / h;
                Geodesic g = c.ray(sx, sy);
                Real L = tracePath(sc.kerr, sc.disk, sc.sky, sc.prop, g, lambda, rng, sc.maxBounces);
                if (L > 0) acc += spec::cieXYZ(lambda).y * L;
            }
            lum.push_back(acc / spp);
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
static void printControls() {
    std::printf(
"\nControls\n"
"  left drag            orbit (inclination and azimuth)\n"
"  right drag           pan the aim point\n"
"  wheel                dolly in/out\n"
"  ctrl + wheel         field of view\n"
"  [ / ]                spin a/M down / up\n"
"  , / .                disk outer radius\n"
"  - / =                exposure\n"
"  b                    cycle scattering bounces (0-4)\n"
"  r                    reset the camera\n"
"  s                    save a PNG snapshot\n"
"  esc or q             quit\n\n");
}

int main(int argc, char** argv) {
    Shared S;
    Scene& sc = S.scene;

    int  winW = 1280, winH = 720;
    Real diskIn = -1, diskOut = 18.0;
    Real targetFps = 30.0;
    double autoQuit = 0.0;      // scripted run: render for N seconds, snapshot, exit
    bool   autoOrbit = false;   // scripted camera motion, to exercise the moving path
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
        else if (a == "--threads")    threads = nextI();
        else if (a == "--autoquit")   autoQuit = nextF();
        else if (a == "--autoorbit")  autoOrbit = true;
        else if (a == "--rtol")       sc.prop.rtol = nextF();
        else if (a == "--help" || a == "-h") {
            std::printf(
"kerrview -- interactive progressive viewer for the Kerr path tracer\n\n"
"  --width N --height N   window size            (1280 720)\n"
"  --spin A               spin a/M               (0.94)\n"
"  --dist R --inc DEG --fov DEG                  (40, 80, 40)\n"
"  --rin R --rout R --tpeak K --albedo A --edge W --turbulence F\n"
"  --bounces N --sky-gain F --nostars\n"
"  --exposure F --fps F --threads N --rtol F\n");
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
    SDL_Renderer* ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
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
    resetAccumulation(S);

    std::printf("kerrview -- %d x %d, %d worker threads\n", winW, winH, nThreads);
    std::printf("  spin a/M %.3f  horizon %.3f M  ISCO %.3f M  disk %.2f -> %.1f M\n",
                sc.kerr.a, sc.kerr.horizon(), sc.kerr.isco(true), sc.disk.rIn, sc.disk.rOut);
    printControls();
    std::fflush(stdout);

    std::vector<std::thread> pool;
    pool.reserve(nThreads);
    for (int i = 0; i < nThreads; ++i) pool.emplace_back(workerLoop, std::ref(S));

    const CameraParams home = S.cam;

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
    bool  running = true;
    int   snapshotIndex = 0;

    auto startTime = Clock::now();
    while (running) {
        // Meter first, on whatever the workers produced since the last reset.
        // Doing this after the reset below always reads an empty buffer, which
        // is exactly how the exposure used to end up meaningless while moving.
        logExposure = updateExposure(S, key, userExposure, logExposure);

        bool camChanged = false;

        // Scripted motion, so the interactive path can be exercised headlessly.
        if (autoOrbit && since(startTime) < autoQuit * 0.6) {
            S.cam.phiDeg += 1.5;
            camChanged = true;
        }
        if (autoQuit > 0 && since(startTime) > autoQuit) {
            std::vector<uint8_t> copy;
            reconfigure(S, [&] { copy = S.display; });
            img::writePng("kerrview-auto.png", winW, winH, copy);
            std::printf("  autoquit: wrote kerrview-auto.png at %u spp, scale 1/%d, %.2f Mrays/s\n",
                        S.minCount.load(), S.rt.scale, raysPerSec / 1e6);
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
                if (dragOrbit && (e.motion.xrel || e.motion.yrel)) {
                    S.cam.phiDeg += e.motion.xrel * 0.25;
                    S.cam.incDeg = clampf(S.cam.incDeg - e.motion.yrel * 0.25, 1.0, 179.0);
                    camChanged = true;
                } else if (dragPan && (e.motion.xrel || e.motion.yrel)) {
                    S.cam.yawDeg   = clampf(S.cam.yawDeg   - e.motion.xrel * 0.05, -80.0, 80.0);
                    S.cam.pitchDeg = clampf(S.cam.pitchDeg + e.motion.yrel * 0.05, -80.0, 80.0);
                    camChanged = true;
                }
                break;

            case SDL_MOUSEWHEEL: {
                bool ctrl = (SDL_GetModState() & KMOD_CTRL) != 0;
                if (ctrl) {
                    S.cam.fovDeg = clampf(S.cam.fovDeg * std::exp(-e.wheel.y * 0.08), 2.0, 120.0);
                } else {
                    Real minR = sc.kerr.horizon() * 2.0 + 1.0;
                    S.cam.camR = clampf(S.cam.camR * std::exp(-e.wheel.y * 0.10), minR, 4000.0);
                    sc.prop.rEscape = std::max<Real>(2000.0, S.cam.camR * 8);
                }
                camChanged = true;
                break;
            }

            case SDL_KEYDOWN:
                switch (e.key.keysym.sym) {
                case SDLK_ESCAPE: case SDLK_q: running = false; break;
                case SDLK_r: S.cam = home; camChanged = true; break;
                case SDLK_LEFTBRACKET:
                case SDLK_RIGHTBRACKET: {
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
                    char name[128];
                    std::snprintf(name, sizeof name, "kerrview-%03d.png", snapshotIndex++);
                    std::vector<uint8_t> copy;
                    reconfigure(S, [&] { copy = S.display; });
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

        // --- choose the regime -------------------------------------------
        bool moving = since(lastInput) < 0.25;
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

        if (camChanged || desiredScale != S.rt.scale) {
            reconfigure(S, [&] {
                S.rt.scale = desiredScale;
                resetAccumulation(S);
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
        SDL_UpdateTexture(tex, nullptr, S.display.data(), winW * 3);
        SDL_RenderCopy(ren, tex, nullptr, nullptr);
        SDL_RenderPresent(ren);

        double dtFrame = since(lastFrame);
        lastFrame = Clock::now();
        if (dtFrame > 0) fps += (1.0 / dtFrame - fps) * 0.1;

        if (autoQuit > 0 && since(lastStatus) > 0.5) {
            double passSamples = double(S.rt.rw) * S.rt.rh;
            double passHz = raysPerSec / std::max(1.0, passSamples);
            std::printf("  t=%4.1fs  %-8s scale 1/%-2d  render %dx%d  %2u spp  ui %.0f fps  "
                        "first-pass %.1f Hz  %.2f Mrays/s\n",
                        since(startTime), moving ? "moving" : "settled", S.rt.scale,
                        S.rt.rw, S.rt.rh, S.minCount.load(), fps, passHz, raysPerSec / 1e6);
            std::printf("           metering: covered %5.1f%%  anchor %.4g  target %.4g  live %.4g\n",
                        100.0 * S.meterCovered, S.meterAnchor, S.meterTarget, S.exposure.load());
            std::fflush(stdout);
            lastStatus = Clock::now();
        }
        if (since(lastTitle) > 0.2) {
            char title[256];
            std::snprintf(title, sizeof title,
                "kerrview  |  %s  1/%d res  |  %u spp  |  %.0f fps  |  %.2f Mrays/s  |  "
                "a=%.2f  r=%.0fM  inc=%.0f  fov=%.0f",
                moving ? "moving" : "settling", S.rt.scale, S.minCount.load(),
                fps, raysPerSec / 1e6,
                sc.kerr.a, double(S.cam.camR), double(S.cam.incDeg), double(S.cam.fovDeg));
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

    SDL_DestroyTexture(tex);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}
