#!/usr/bin/env bash
# Build the Kerr path tracer and the interactive viewer on Linux.
#
# The Unix counterpart of build.ps1, and deliberately the same shape:
#
#   ./build.sh               both targets
#   ./build.sh cli           batch renderer only  -> kerr
#   ./build.sh viewer        interactive viewer   -> kerrview
#
# The source is portable -- neither front end has a line of platform-specific
# code in it -- so all this script does is find a compiler and, for the viewer,
# find SDL2.
#
# Compilers are tried in order: $CXX, clang++, g++, then a copy of Zig under
# .toolchain, fetched on first use.  `zig c++` is a self-contained Clang and
# needs no admin install, which is what makes the fallback worth having; if a
# system compiler is already present it is used instead and nothing is
# downloaded.
#
# SDL2 comes from the distribution, since unlike Windows there is a packaged
# one:
#
#   sudo apt install libsdl2-dev        # Debian / Ubuntu
#
# Any recent Clang or GCC works without this script at all:
#   clang++ -std=c++20 -O3 -march=native src/main.cpp -o kerr

set -euo pipefail

target="${1:-all}"
case "$target" in
    all|cli|viewer) ;;
    -h|--help) sed -n '2,27p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown target '$target' (expected: all, cli, viewer)" >&2; exit 2 ;;
esac

root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
toolDir="$root/.toolchain"
zigVer='0.16.0'

case "$(uname -m)" in
    x86_64|amd64)  zigArch='x86_64'  ;;
    aarch64|arm64) zigArch='aarch64' ;;
    *)             zigArch="$(uname -m)" ;;
esac
zigDir="$toolDir/zig-${zigArch}-linux-${zigVer}"

# Compiler probes need care.  The source has to be a real file, because
# `zig c++` will not read /dev/null, and the probe has to compile to a real
# object, because `zig c++ -fsyntax-only` is broken in Zig 0.16 -- it fails with
# "FileNotFound" whatever it is handed.  Either mistake makes every probe below
# report failure and the build quietly comes out wrong.
probe="$(mktemp "${TMPDIR:-/tmp}/kerr-probe.XXXXXX.cpp")"
probeObj="${probe%.cpp}.o"
trap 'rm -f "$probe" "$probeObj"' EXIT
printf 'int main() { return 0; }\n' > "$probe"

# ---------------------------------------------------------------------------
# Compiler
# ---------------------------------------------------------------------------
# Read $CXX before shadowing it.  The compiler is held as an array because the
# Zig fallback is two words, `zig c++`.
cxxEnv="${CXX:-}"
declare -a CXX=()

works() { "$@" -std=c++20 -c "$probe" -o "$probeObj" >/dev/null 2>&1; }

get_zig() {
    if [ -x "$zigDir/zig" ]; then echo "$zigDir/zig"; return; fi
    if command -v zig >/dev/null 2>&1; then command -v zig; return; fi

    local url="https://ziglang.org/download/${zigVer}/zig-${zigArch}-linux-${zigVer}.tar.xz"
    echo "fetching Zig toolchain (~55 MB) ..." >&2
    mkdir -p "$toolDir"
    local tarball="$toolDir/zig.tar.xz"
    if command -v curl >/dev/null 2>&1; then
        curl -fsSL -o "$tarball" "$url"
    elif command -v wget >/dev/null 2>&1; then
        wget -q -O "$tarball" "$url"
    else
        echo "no C++ compiler found, and neither curl nor wget is available to fetch one" >&2
        exit 1
    fi
    tar -xJf "$tarball" -C "$toolDir"
    rm -f "$tarball"
    echo "$zigDir/zig"
}

if [ -n "$cxxEnv" ]; then
    read -r -a CXX <<< "$cxxEnv"
    works "${CXX[@]}" || { echo "CXX='$cxxEnv' cannot compile C++20" >&2; exit 1; }
elif command -v clang++ >/dev/null 2>&1 && works clang++; then
    CXX=(clang++)
elif command -v g++ >/dev/null 2>&1 && works g++; then
    CXX=(g++)
else
    CXX=("$(get_zig)" c++)
    works "${CXX[@]}" || { echo "the Zig toolchain in .toolchain cannot compile C++20" >&2; exit 1; }
fi
echo "compiler: ${CXX[*]}"

# ---------------------------------------------------------------------------
# Flags
# ---------------------------------------------------------------------------
flags=(
    -std=c++20 -O3
    -funroll-loops
    -fno-math-errno
    -fno-trapping-math
    -ffp-contract=fast
    -fomit-frame-pointer
    -pthread
    -Wno-nullability-completeness
)
# Deliberately NOT -ffast-math: the integrator uses isfinite() to reject bad
# trial steps, and -ffinite-math-only would compile those checks away.

# -march=native is worth ~16%, and it is FMA that earns it, not vector width --
# the hot loop is scalar doubles.  It is probed rather than assumed, because it
# is not accepted by every compiler on every architecture (notably older
# aarch64 GCC), and a build that works beats the last 16%.
if works "${CXX[@]}" -march=native; then
    flags+=(-march=native)
else
    echo "note: -march=native is not accepted here; building without it"
fi

# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------
# Libraries go after the translation unit, not before it: GNU ld links with
# --as-needed by default on Debian and Ubuntu, and drops a -l that appears
# before the object that references it.
declare -a extraCFlags=() extraLibs=()

build() {
    local name="$1"
    local out; [ "$name" = 'main' ] && out="$root/kerr" || out="$root/kerrview"

    echo "compiling -> $out"
    local start end; start=$(date +%s%N)
    "${CXX[@]}" "${flags[@]}" ${extraCFlags[@]+"${extraCFlags[@]}"} \
        "$root/src/$name.cpp" -o "$out" ${extraLibs[@]+"${extraLibs[@]}"}
    end=$(date +%s%N)
    LC_ALL=C awk -v ns="$((end - start))" 'BEGIN { printf "  done in %.1fs\n", ns / 1e9 }'
}

if [ "$target" = 'all' ] || [ "$target" = 'cli' ]; then
    extraCFlags=(); extraLibs=()
    build main
fi

if [ "$target" = 'all' ] || [ "$target" = 'viewer' ]; then
    if command -v pkg-config >/dev/null 2>&1 && pkg-config --exists sdl2; then
        # pkg-config rather than a hardcoded -I/usr/include/SDL2: Debian and
        # Ubuntu split the headers, keeping SDL_config.h in the multiarch
        # include directory, and only pkg-config knows about both halves.
        read -r -a extraCFlags <<< "$(pkg-config --cflags sdl2)"
        read -r -a extraLibs   <<< "$(pkg-config --libs sdl2)"
    elif command -v sdl2-config >/dev/null 2>&1; then
        read -r -a extraCFlags <<< "$(sdl2-config --cflags)"
        read -r -a extraLibs   <<< "$(sdl2-config --libs)"
        # sdl2-config offers only <prefix>/include/SDL2, but the viewer includes
        # <SDL2/SDL.h>, so <prefix>/include has to be on the path as well.  That
        # is free on a system install and necessary for any other prefix.
        extraCFlags+=(-I"$(sdl2-config --prefix)/include")
    else
        echo "SDL2 development files not found -- the viewer needs them." >&2
        echo "  Debian/Ubuntu:  sudo apt install libsdl2-dev" >&2
        echo "  Fedora:         sudo dnf install SDL2-devel" >&2
        echo "  Arch:           sudo pacman -S sdl2" >&2
        echo "(./build.sh cli builds the batch renderer, which needs nothing.)" >&2
        exit 1
    fi
    # Debian and Ubuntu split the headers: the generated SDL_config.h lives in
    # the multiarch include directory and the plain one only forwards to it, so
    # <SDL2/_real_SDL_config.h> has to be reachable or nothing compiles.
    # Neither pkg-config nor sdl2-config names that directory.  A compiler that
    # knows the distribution adds it on its own -- `zig c++` does, once it
    # exists -- so a packaged SDL2 needs nothing from this loop; an SDL2 under
    # any other prefix does.  Look beside each include directory we were handed
    # rather than assuming /usr, which is the case that is actually at risk.
    declare -a multiarch=()
    for f in ${extraCFlags[@]+"${extraCFlags[@]}"}; do
        [ "${f#-I}" = "$f" ] && continue
        for d in "${f#-I}"/*-linux-*; do
            if [ -d "$d/SDL2" ]; then multiarch+=(-I"$d"); fi
        done
    done
    extraCFlags+=(${multiarch[@]+"${multiarch[@]}"})

    build viewer
fi

echo ""
echo "run: ./kerr --check      validate the geodesic engine"
echo "     ./kerr --preview    fast low-res still"
echo "     ./kerrview          interactive viewer"
