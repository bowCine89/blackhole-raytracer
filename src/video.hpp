// video.hpp -- streaming uncompressed AVI writer.
//
// The rest of the project has no dependencies and this keeps it that way.  An
// uncompressed RIFF/AVI holds the frames exactly as rendered -- no colour
// subsampling, no block artefacts, nothing between the tone mapper and the file
// -- which is the point when the frames each took seconds to converge.  The
// cost is size: 1280x720x3 bytes a frame, so about 265 MB for a 96-frame loop.
// Transcode with ffmpeg afterwards if something smaller is wanted; going the
// other way is impossible.
//
// Frames are written as they arrive rather than buffered, so the batch renderer
// never holds the whole sequence in memory.  The header counts are patched on
// close, once the frame count is known.
#pragma once
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#if defined(_WIN32)
#  define KERR_POPEN   _popen
#  define KERR_PCLOSE  _pclose
#  define KERR_PIPEW   "wb"
#  define KERR_NULLDEV "NUL"
#else
#  define KERR_POPEN   popen
#  define KERR_PCLOSE  pclose
#  define KERR_PIPEW   "w"
#  define KERR_NULLDEV "/dev/null"
#  include <csignal>
#endif

namespace vid {

class AviWriter {
public:
    bool open(const std::string& path, int width, int height, double fps) {
        w_ = width; h_ = height; fps_ = fps > 0 ? fps : 24.0;
        stride_ = ((w_ * 3 + 3) / 4) * 4;          // AVI rows are 4-byte aligned
        frameBytes_ = size_t(stride_) * h_;
        row_.resize(stride_, 0);
        f_ = std::fopen(path.c_str(), "wb");
        if (!f_) return false;
        writeHeaders(0);                            // frame count patched in close()
        moviListPos_ = tell();
        put("LIST"); u32(0); put("movi");           // size patched in close()
        moviDataPos_ = moviListPos_ + 8;            // idx1 offsets are relative to 'movi'
        return true;
    }

    // rgb is width*height*3, row-major, top row first -- what the renderer has.
    // AVI wants BGR, bottom row first, so both are flipped on the way out.
    bool addFrame(const uint8_t* rgb) {
        if (!f_) return false;
        uint32_t chunkPos = tell();
        put("00db");
        u32(uint32_t(frameBytes_));
        for (int y = h_ - 1; y >= 0; --y) {
            const uint8_t* src = rgb + size_t(y) * w_ * 3;
            uint8_t* dst = row_.data();
            for (int x = 0; x < w_; ++x) {
                dst[0] = src[2]; dst[1] = src[1]; dst[2] = src[0];   // RGB -> BGR
                dst += 3; src += 3;
            }
            std::fwrite(row_.data(), 1, stride_, f_);
        }
        if (frameBytes_ & 1) u8(0);                 // RIFF chunks are word aligned
        index_.push_back(chunkPos - moviDataPos_);
        ++frames_;
        return true;
    }

    bool close() {
        if (!f_) return false;
        uint32_t moviEnd = tell();

        // idx1: one entry per frame, each a key frame.
        put("idx1");
        u32(uint32_t(index_.size() * 16));
        for (uint32_t off : index_) {
            put("00db");
            u32(0x10);                              // AVIIF_KEYFRAME
            u32(off);
            u32(uint32_t(frameBytes_));
        }
        uint32_t fileEnd = tell();

        // Patch the sizes now that they are known.  The header rewrite has to
        // come first: it starts with a fresh `RIFF` + placeholder size, so
        // patching offset 4 before it would just be overwritten with zero.
        seek(0);                writeHeaders(frames_);              // header frame counts
        seek(4);                u32(fileEnd - 8);                   // RIFF size
        seek(moviListPos_ + 4); u32(moviEnd - moviListPos_ - 8);    // movi LIST size

        std::fclose(f_);
        f_ = nullptr;
        return true;
    }

    int frames() const { return frames_; }
    size_t bytes() const { return size_t(frames_) * frameBytes_; }

private:
    void put(const char* four) { std::fwrite(four, 1, 4, f_); }
    void u8(uint8_t v)  { std::fwrite(&v, 1, 1, f_); }
    void u16(uint16_t v) { uint8_t b[2] = {uint8_t(v), uint8_t(v >> 8)}; std::fwrite(b, 1, 2, f_); }
    void u32(uint32_t v) {
        uint8_t b[4] = {uint8_t(v), uint8_t(v >> 8), uint8_t(v >> 16), uint8_t(v >> 24)};
        std::fwrite(b, 1, 4, f_);
    }
    uint32_t tell() { return uint32_t(std::ftell(f_)); }
    void seek(uint32_t p) { std::fseek(f_, long(p), SEEK_SET); }

    void writeHeaders(int frameCount) {
        const uint32_t usPerFrame = uint32_t(1.0e6 / fps_ + 0.5);

        put("RIFF"); u32(0); put("AVI ");

        put("LIST"); u32(4 + 8 + 56 + 4 + 8 + 8 + 56 + 8 + 40); put("hdrl");

        put("avih"); u32(56);
        u32(usPerFrame);
        u32(uint32_t(frameBytes_ * fps_));      // max bytes/sec
        u32(0);                                  // padding granularity
        u32(0x10);                               // AVIF_HASINDEX
        u32(uint32_t(frameCount));
        u32(0);                                  // initial frames
        u32(1);                                  // streams
        u32(uint32_t(frameBytes_));              // suggested buffer
        u32(uint32_t(w_)); u32(uint32_t(h_));
        u32(0); u32(0); u32(0); u32(0);          // reserved

        put("LIST"); u32(4 + 8 + 56 + 8 + 40); put("strl");

        put("strh"); u32(56);
        put("vids"); put("DIB ");                // uncompressed device-independent bitmap
        u32(0); u16(0); u16(0);
        u32(0);                                  // initial frames
        u32(1000);                               // scale
        u32(uint32_t(fps_ * 1000.0 + 0.5));      // rate: fractional fps stays exact
        u32(0);                                  // start
        u32(uint32_t(frameCount));               // length
        u32(uint32_t(frameBytes_));
        u32(0xFFFFFFFFu);                        // quality: default
        u32(0);                                  // sample size (0 = variable)
        u16(0); u16(0); u16(uint16_t(w_)); u16(uint16_t(h_));

        put("strf"); u32(40);
        u32(40);                                 // BITMAPINFOHEADER size
        u32(uint32_t(w_));
        u32(uint32_t(h_));                       // positive: rows bottom-up
        u16(1); u16(24);                         // planes, bits per pixel
        u32(0);                                  // BI_RGB, no compression
        u32(uint32_t(frameBytes_));
        u32(0); u32(0); u32(0); u32(0);
    }

    FILE* f_ = nullptr;
    int w_ = 0, h_ = 0, stride_ = 0, frames_ = 0;
    double fps_ = 24.0;
    size_t frameBytes_ = 0;
    uint32_t moviListPos_ = 0, moviDataPos_ = 0;
    std::vector<uint32_t> index_;
    std::vector<uint8_t> row_;
};

// ---------------------------------------------------------------------------
// H.265, by piping raw frames to ffmpeg
// ---------------------------------------------------------------------------
//
// An HEVC encoder is not something that belongs in a header, and linking libx265
// would end the project's habit of building with nothing installed.  So ffmpeg
// is a *runtime* dependency and an optional one: frames go down a pipe as raw
// RGB and it does the encoding.  The build is unchanged, and a machine without
// ffmpeg still renders everything -- it just cannot compress it.
//
// The defaults are chosen for this material rather than for video in general.
// The disk is a smooth gradient over a black field, which is exactly what shows
// 8-bit banding, so the pipe is 10-bit by default even though the source is
// 8-bit RGB: the extra depth costs almost nothing and gives the encoder room to
// dither rather than contour.  crf 18 at preset slow is visually lossless on
// this content; crf 0 is mathematically lossless if that is what you want.
class H265Writer {
public:
    int         crf    = 18;
    std::string preset = "slow";
    bool        tenBit = true;

    // ffmpeg is looked up on PATH unless KERR_FFMPEG names a binary.
    static std::string exe() {
        const char* e = std::getenv("KERR_FFMPEG");
        return (e && *e) ? std::string(e) : std::string("ffmpeg");
    }

    static bool available() {
        std::string cmd = quote(exe()) + " -hide_banner -version > " KERR_NULLDEV " 2>&1";
        return std::system(cmd.c_str()) == 0;
    }

    bool open(const std::string& path, int width, int height, double fps) {
        w_ = width; h_ = height; fps_ = fps > 0 ? fps : 24.0;
        path_ = path;
        frameBytes_ = size_t(w_) * h_ * 3;

        char buf[1024];
        std::snprintf(buf, sizeof buf,
            "%s -hide_banner -loglevel error -y -f rawvideo -pix_fmt rgb24 "
            "-s %dx%d -r %.6f -i - -an -c:v libx265 -preset %s -crf %d "
            "-pix_fmt %s -x265-params log-level=error -tag:v hvc1 %s",
            quote(exe()).c_str(), w_, h_, fps_, preset.c_str(), crf,
            tenBit ? "yuv420p10le" : "yuv420p", quote(path_).c_str());

#if !defined(_WIN32)
        // If the encoder dies mid-sequence the next write hits a broken pipe,
        // and the default disposition for SIGPIPE would take the renderer down
        // with it -- silently, halfway through a job that can run for hours.
        // Ignoring it turns that into a short write we can report instead.
        std::signal(SIGPIPE, SIG_IGN);
#endif
        f_ = KERR_POPEN(buf, KERR_PIPEW);
        return f_ != nullptr;
    }

    bool addFrame(const uint8_t* rgb) {
        if (!f_) return false;
        if (std::fwrite(rgb, 1, frameBytes_, f_) != frameBytes_) {
            broken_ = true;          // encoder went away; close() will report it
            return false;
        }
        ++frames_;
        return true;
    }

    bool close() {
        if (!f_) return false;
        int rc = KERR_PCLOSE(f_);
        f_ = nullptr;
        // Size comes off the finished file: the whole point is that we do not
        // know it in advance.
        if (FILE* g = std::fopen(path_.c_str(), "rb")) {
            std::fseek(g, 0, SEEK_END);
            long n = std::ftell(g);
            std::fclose(g);
            if (n > 0) bytes_ = size_t(n);
        }
        return rc == 0 && !broken_;
    }

    int    frames() const { return frames_; }
    size_t bytes()  const { return bytes_; }

private:
    // Shell-quote a path.  Single quotes everywhere but Windows, which does not
    // understand them.
    static std::string quote(const std::string& s) {
#if defined(_WIN32)
        return "\"" + s + "\"";
#else
        std::string o = "'";
        for (char c : s) { if (c == '\'') o += "'\\''"; else o += c; }
        return o + "'";
#endif
    }

    FILE* f_ = nullptr;
    std::string path_;
    int w_ = 0, h_ = 0, frames_ = 0;
    double fps_ = 24.0;
    size_t frameBytes_ = 0, bytes_ = 0;
    bool   broken_ = false;
};

// ---------------------------------------------------------------------------
// One handle over both, chosen by the output extension
// ---------------------------------------------------------------------------
// .avi keeps the built-in uncompressed writer, so nothing that worked before
// behaves differently.  Anything ffmpeg can mux -- .mp4, .mkv, .mov, .hevc --
// goes through H.265.
class VideoWriter {
public:
    int crf = 18;

    static bool wantsH265(const std::string& path) {
        size_t dot = path.find_last_of('.');
        if (dot == std::string::npos) return false;
        std::string e = path.substr(dot + 1);
        for (char& c : e) c = char(std::tolower((unsigned char)c));
        return e == "mp4" || e == "mkv" || e == "mov" || e == "hevc" || e == "265";
    }

    // Returns false with a printed reason; the caller decides what to do.
    bool open(const std::string& path, int w, int h, double fps) {
        h265_ = wantsH265(path);
        if (!h265_) return avi_.open(path, w, h, fps);
        if (!H265Writer::available()) {
            std::fprintf(stderr,
                "  %s needs ffmpeg with libx265 on PATH, and it was not found.\n"
                "    Debian/Ubuntu:  sudo apt install ffmpeg\n"
                "    or set KERR_FFMPEG to a binary, or write .avi instead.\n",
                path.c_str());
            return false;
        }
        hev_.crf = crf;
        return hev_.open(path, w, h, fps);
    }
    bool addFrame(const uint8_t* rgb) { return h265_ ? hev_.addFrame(rgb) : avi_.addFrame(rgb); }
    bool close()                      { return h265_ ? hev_.close()       : avi_.close(); }
    int    frames() const             { return h265_ ? hev_.frames()      : avi_.frames(); }
    size_t bytes()  const             { return h265_ ? hev_.bytes()       : avi_.bytes(); }
    bool   compressed() const         { return h265_; }
    const char* codec() const         { return h265_ ? "H.265" : "uncompressed"; }

private:
    bool h265_ = false;
    AviWriter  avi_;
    H265Writer hev_;
};

} // namespace vid
