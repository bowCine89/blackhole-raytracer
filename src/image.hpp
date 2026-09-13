// image.hpp -- dependency-free image output (PNG, PPM, PFM).
//
// The PNG writer emits a valid zlib stream built from uncompressed DEFLATE
// blocks.  That costs file size but keeps the whole project free of external
// libraries, which matters more here than a few megabytes on disk.
#pragma once
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace img {

inline uint32_t crcTable(int n) {
    static uint32_t table[256];
    static bool built = false;
    if (!built) {
        for (uint32_t i = 0; i < 256; ++i) {
            uint32_t c = i;
            for (int k = 0; k < 8; ++k) c = (c & 1) ? 0xEDB88320u ^ (c >> 1) : (c >> 1);
            table[i] = c;
        }
        built = true;
    }
    return table[n];
}

inline uint32_t crc32Buf(const uint8_t* d, size_t n, uint32_t crc = 0xFFFFFFFFu) {
    for (size_t i = 0; i < n; ++i) crc = crcTable((crc ^ d[i]) & 0xFF) ^ (crc >> 8);
    return crc;
}

inline uint32_t adler32Buf(const uint8_t* d, size_t n) {
    uint32_t a = 1, b = 0;
    for (size_t i = 0; i < n; ++i) { a = (a + d[i]) % 65521; b = (b + a) % 65521; }
    return (b << 16) | a;
}

inline void put32(std::vector<uint8_t>& v, uint32_t x) {
    v.push_back(uint8_t(x >> 24)); v.push_back(uint8_t(x >> 16));
    v.push_back(uint8_t(x >> 8));  v.push_back(uint8_t(x));
}

inline void chunk(std::vector<uint8_t>& out, const char* type, const std::vector<uint8_t>& data) {
    put32(out, uint32_t(data.size()));
    size_t start = out.size();
    out.insert(out.end(), type, type + 4);
    out.insert(out.end(), data.begin(), data.end());
    uint32_t c = crc32Buf(out.data() + start, out.size() - start) ^ 0xFFFFFFFFu;
    put32(out, c);
}

// rgb8 is width*height*3 bytes, row-major, top row first.
inline bool writePng(const std::string& path, int w, int h, const std::vector<uint8_t>& rgb8) {
    std::vector<uint8_t> raw;
    raw.reserve(size_t(h) * (size_t(w) * 3 + 1));
    for (int y = 0; y < h; ++y) {
        raw.push_back(0);                               // filter type: none
        const uint8_t* row = rgb8.data() + size_t(y) * w * 3;
        raw.insert(raw.end(), row, row + size_t(w) * 3);
    }

    std::vector<uint8_t> z;
    z.push_back(0x78); z.push_back(0x01);               // zlib header, no preset dict
    size_t pos = 0;
    while (pos < raw.size()) {
        size_t n = std::min<size_t>(65535, raw.size() - pos);
        bool last = (pos + n == raw.size());
        z.push_back(last ? 1 : 0);                      // stored block
        z.push_back(uint8_t(n & 0xFF));       z.push_back(uint8_t(n >> 8));
        z.push_back(uint8_t(~n & 0xFF));      z.push_back(uint8_t((~n >> 8) & 0xFF));
        z.insert(z.end(), raw.begin() + pos, raw.begin() + pos + n);
        pos += n;
    }
    put32(z, adler32Buf(raw.data(), raw.size()));

    std::vector<uint8_t> out = {137, 80, 78, 71, 13, 10, 26, 10};
    std::vector<uint8_t> ihdr;
    put32(ihdr, uint32_t(w)); put32(ihdr, uint32_t(h));
    ihdr.push_back(8); ihdr.push_back(2); ihdr.push_back(0); ihdr.push_back(0); ihdr.push_back(0);
    chunk(out, "IHDR", ihdr);
    chunk(out, "IDAT", z);
    chunk(out, "IEND", {});

    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    bool ok = std::fwrite(out.data(), 1, out.size(), f) == out.size();
    std::fclose(f);
    return ok;
}

inline bool writePpm(const std::string& path, int w, int h, const std::vector<uint8_t>& rgb8) {
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    std::fprintf(f, "P6\n%d %d\n255\n", w, h);
    bool ok = std::fwrite(rgb8.data(), 1, rgb8.size(), f) == rgb8.size();
    std::fclose(f);
    return ok;
}

// 32-bit float RGB, little-endian, bottom row first (PFM convention).
inline bool writePfm(const std::string& path, int w, int h, const std::vector<float>& rgb) {
    FILE* f = std::fopen(path.c_str(), "wb");
    if (!f) return false;
    std::fprintf(f, "PF\n%d %d\n-1.0\n", w, h);
    bool ok = true;
    for (int y = h - 1; y >= 0 && ok; --y)
        ok = std::fwrite(rgb.data() + size_t(y) * w * 3, sizeof(float), size_t(w) * 3, f)
             == size_t(w) * 3;
    std::fclose(f);
    return ok;
}

} // namespace img
