# Build the Kerr path tracer.
#
# Uses the portable Zig toolchain in .toolchain (zig c++ is Clang + libc++,
# self-contained, no installer and no admin rights).  Any recent clang++ or
# g++ works just as well:
#
#   clang++ -std=c++20 -O3 -march=native src/main.cpp -o kerr.exe
#   g++     -std=c++20 -O3 -march=native src/main.cpp -o kerr

$zig = Join-Path $PSScriptRoot '.toolchain\zig-x86_64-windows-0.16.0\zig.exe'
if (-not (Test-Path $zig)) {
    $cmd = Get-Command zig -ErrorAction SilentlyContinue
    if ($cmd) { $zig = $cmd.Source }
    else {
        Write-Error "No toolchain found. Expected .toolchain\zig-x86_64-windows-0.16.0\zig.exe, or zig / clang++ on PATH."
        exit 1
    }
}

$flags = @(
    '-std=c++20', '-O3',
    '-march=native',          # worth ~16%, and it is FMA that earns it, not
                              # vector width -- the hot loop is scalar doubles
    '-funroll-loops',
    '-fno-math-errno',        # sin/cos/sqrt need not set errno
    '-fno-trapping-math',
    '-ffp-contract=fast',     # allow FMA contraction
    '-fomit-frame-pointer',
    '-Wno-nullability-completeness'
)
# Deliberately NOT -ffast-math: the integrator uses isfinite() to reject bad
# trial steps, and -ffinite-math-only would compile those checks away.

$src = Join-Path $PSScriptRoot 'src\main.cpp'
$out = Join-Path $PSScriptRoot 'kerr.exe'

Write-Host "compiling -> $out"
$sw = [Diagnostics.Stopwatch]::StartNew()

# Link to a scratch name and copy into place.  With Smart App Control in
# enforcement mode, Windows sometimes blocks an unsigned binary that a compiler
# has just written in place ("an application control policy has blocked this
# file"); writing fresh and copying avoids that verdict.
$tmp = Join-Path $PSScriptRoot 'kerr.build.exe'

# Native tools write progress to stderr; don't let PowerShell read that as a
# failure.  The exit code is the authority.
$prev = $ErrorActionPreference
$ErrorActionPreference = 'Continue'
& $zig c++ @flags $src -o $tmp 2>&1 | ForEach-Object { "$_" } | Write-Host
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
Remove-Item (Join-Path $PSScriptRoot 'kerr.build.pdb') -ErrorAction SilentlyContinue
Write-Host ("done in {0:N1}s" -f $sw.Elapsed.TotalSeconds)
Write-Host "run: .\kerr.exe --check      (validate the geodesic engine)"
Write-Host "     .\kerr.exe --preview    (fast low-res image)"
Write-Host "     .\kerr.exe              (1280x720, 128 spp)"
