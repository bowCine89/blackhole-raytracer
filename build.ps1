# Build the Kerr path tracer and the interactive viewer.
#
# Both targets come from one source tree: src/render.hpp holds the tracer, and
# main.cpp / viewer.cpp are two front ends onto it.
#
#   .\build.ps1              both targets
#   .\build.ps1 cli          batch renderer only  -> kerr.exe
#   .\build.ps1 viewer       interactive viewer   -> kerrview.exe
#
# Dependencies live under .toolchain and are fetched on first use:
#   * Zig     -- `zig c++` is a self-contained Clang + libc++, no admin install
#   * SDL2    -- MinGW development package, for the viewer only
#
# Any recent Clang or GCC works instead of Zig:
#   clang++ -std=c++20 -O3 -march=native src/main.cpp -o kerr

param([string]$Target = "all")

$ErrorActionPreference = 'Stop'
$root    = $PSScriptRoot
$toolDir = Join-Path $root '.toolchain'
$zigDir  = Join-Path $toolDir 'zig-x86_64-windows-0.16.0'
$sdlVer  = '2.30.9'
$sdlDir  = Join-Path $toolDir "SDL2-$sdlVer\x86_64-w64-mingw32"

function Get-Zig {
    $exe = Join-Path $zigDir 'zig.exe'
    if (Test-Path $exe) { return $exe }
    $cmd = Get-Command zig -ErrorAction SilentlyContinue
    if ($cmd) { return $cmd.Source }
    Write-Host "fetching Zig toolchain (~93 MB) ..."
    New-Item -ItemType Directory -Force $toolDir | Out-Null
    $zip = Join-Path $toolDir 'zig.zip'
    (New-Object Net.WebClient).DownloadFile(
        "https://ziglang.org/download/0.16.0/zig-x86_64-windows-0.16.0.zip", $zip)
    Expand-Archive -Path $zip -DestinationPath $toolDir -Force
    Remove-Item $zip
    return (Join-Path $zigDir 'zig.exe')
}

function Get-Sdl {
    if (Test-Path (Join-Path $sdlDir 'include\SDL2\SDL.h')) { return }
    Write-Host "fetching SDL2 $sdlVer (~7 MB) ..."
    New-Item -ItemType Directory -Force $toolDir | Out-Null
    $zip = Join-Path $toolDir 'sdl.zip'
    (New-Object Net.WebClient).DownloadFile(
        "https://github.com/libsdl-org/SDL/releases/download/release-$sdlVer/SDL2-devel-$sdlVer-mingw.zip", $zip)
    Expand-Archive -Path $zip -DestinationPath $toolDir -Force
    Remove-Item $zip
}

$zig = Get-Zig

$flags = @(
    '-std=c++20', '-O3',
    '-march=native',          # worth ~16%, and it is FMA that earns it, not
                              # vector width -- the hot loop is scalar doubles
    '-funroll-loops',
    '-fno-math-errno',
    '-fno-trapping-math',
    '-ffp-contract=fast',
    '-fomit-frame-pointer',
    '-Wno-nullability-completeness'
)
# Deliberately NOT -ffast-math: the integrator uses isfinite() to reject bad
# trial steps, and -ffinite-math-only would compile those checks away.

# Link to a scratch name and copy into place.  With Smart App Control in
# enforcement mode, Windows sometimes blocks an unsigned binary a compiler has
# just written in place; writing fresh and copying avoids that verdict.
function Invoke-Build([string]$name, [string[]]$extra) {
    $src = Join-Path $root "src\$name.cpp"
    $out = Join-Path $root $(if ($name -eq 'main') { 'kerr.exe' } else { 'kerrview.exe' })
    $tmp = Join-Path $root "$name.build.exe"

    Write-Host "compiling -> $out"
    $sw = [Diagnostics.Stopwatch]::StartNew()
    $prev = $ErrorActionPreference
    $ErrorActionPreference = 'Continue'
    & $zig c++ @flags @extra $src -o $tmp 2>&1 | ForEach-Object { "$_" } | Write-Host
    $code = $LASTEXITCODE
    $ErrorActionPreference = $prev
    $sw.Stop()

    if ($code -ne 0) {
        Remove-Item $tmp -ErrorAction SilentlyContinue
        Write-Error "build failed (exit $code)"
        exit $code
    }
    Copy-Item $tmp $out -Force
    Remove-Item $tmp -ErrorAction SilentlyContinue
    Remove-Item (Join-Path $root "$name.build.pdb") -ErrorAction SilentlyContinue
    Write-Host ("  done in {0:N1}s" -f $sw.Elapsed.TotalSeconds)
}

if ($Target -eq 'all' -or $Target -eq 'cli') {
    Invoke-Build 'main' @()
}

if ($Target -eq 'all' -or $Target -eq 'viewer') {
    Get-Sdl
    $sdlFlags = @(
        "-I$sdlDir\include",
        "-L$sdlDir\lib",
        '-lSDL2main', '-lSDL2',
        '-lole32', '-loleaut32', '-limm32', '-lversion',
        '-lsetupapi', '-lwinmm', '-lgdi32', '-luser32', '-lshell32', '-lrpcrt4'
    )
    Invoke-Build 'viewer' $sdlFlags
    Copy-Item (Join-Path $sdlDir 'bin\SDL2.dll') $root -Force   # needed at runtime
}

Write-Host ""
Write-Host "run: .\kerr.exe --check      validate the geodesic engine"
Write-Host "     .\kerr.exe --preview    fast low-res still"
Write-Host "     .\kerrview.exe          interactive viewer"
