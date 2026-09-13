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
#include <cstring>
#include <string>
#include <vector>

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

} // namespace vid
