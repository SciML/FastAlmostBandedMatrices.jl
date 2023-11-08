# FastAlmostBandedMatrices.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)

[![CI](https://github.com/avik-pal/FastAlmostBandedMatrices.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/avik-pal/FastAlmostBandedMatrices.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/LuxDL/FastAlmostBandedMatrices.jl/branch/main/graph/badge.svg?)](https://codecov.io/gh/LuxDL/FastAlmostBandedMatrices.jl)
[![Package Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/FastAlmostBandedMatrices)](https://pkgs.genieframework.com?packages=FastAlmostBandedMatrices)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

A fast implementation of Almost Banded Matrices in Julia. Primarily developed for use in
[BoundaryValueDiffEq.jl](https://github.com/SciML/BoundaryValueDiffEq.jl).

[SemiseparableMatrices.jl](https://github.com/JuliaLinearAlgebra/SemiseparableMatrices.jl)
provides an alternate implementation of almost banded matrices, but that had considerable
overhead for our use-case. This implementation borrows a lot of details from that
repository.

## High Level Examples

```julia
using FastAlmostBandedMatrices

m = 2
n = 10

A1 = AlmostBandedMatrix(brand(Float64, n, n, m + 1, m), rand(Float64, m, n))
```

## Benchmarks

```julia
using BenchmarkTools, FastAlmostBandedMatrices, SparseArrays, FillArrays, LinearAlgebra
import SemiseparableMatrices

m = 5
n = 1000

II = zeros(n, m)
II[diagind(II)] .= 1

A1 = AlmostBandedMatrix(brand(Float64, n, n, m + 1, m), rand(Float64, m, n))
A2 = Matrix(A1)
A3 = SemiseparableMatrices.AlmostBandedMatrix(copy(A1.bands),
    SemiseparableMatrices.LowRankMatrix(II, copy(A1.fill)))
A4 = sparse(A2)

@benchmark qr($A1)
# BenchmarkTools.Trial: 6579 samples with 1 evaluation.
#  Range (min … max):  707.260 μs …   5.453 ms  ┊ GC (min … max): 0.00% … 76.64%
#  Time  (median):     731.708 μs               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   757.683 μs ± 246.801 μs  ┊ GC (mean ± σ):  2.62% ±  6.52%

#    ▆▄█▄                                                          
#   ▅████▇▄▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂▂▁▂▂▂▂▂▂▂▂▂▂▂▂▁▂▂▁▂▂▂▁▁▂▂ ▃
#   707 μs           Histogram: frequency by time         1.15 ms <

#  Memory estimate: 320.70 KiB, allocs estimate: 11.

@benchmark qr($A2)
# BenchmarkTools.Trial: 179 samples with 1 evaluation.
#  Range (min … max):  25.307 ms … 35.189 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     27.488 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   27.975 ms ±  2.168 ms  ┊ GC (mean ± σ):  1.55% ± 3.39%

#   ▂█▄▃ ▃  ▂ ▂ ▄▂ ▂  ▂ ▂▂                                       
#   ██████▇███████▇████▆██▆▆▅▃▆▅▅▇▁▅▃▃▅▇▆▃▃▃▃▃▁▃▁▁▃▁▁▃▃▃▃▁▁▁▅▁▃ ▃
#   25.3 ms         Histogram: frequency by time        34.4 ms <

#  Memory estimate: 8.18 MiB, allocs estimate: 6.

@benchmark qr($A3)
# BenchmarkTools.Trial: 900 samples with 1 evaluation.
#  Range (min … max):  3.149 ms … 16.527 ms  ┊ GC (min … max):  0.00% … 69.59%
#  Time  (median):     4.157 ms              ┊ GC (median):     0.00%
#  Time  (mean ± σ):   5.557 ms ±  3.123 ms  ┊ GC (mean ± σ):  20.55% ± 22.17%

#   ▆█▃▁ ▁       ▂▆▄                                    ▁       
#   ██████▇▇█▇██████▇▄▅▅▁▁▁▄▁▄▁▁▁▁▁▁▄▁▁▁▁▁▁▁▁▆▅▇▇▇██▆█▆██▇▇▇▅▆ █
#   3.15 ms      Histogram: log(frequency) by time       14 ms <

#  Memory estimate: 6.32 MiB, allocs estimate: 65108.

@benchmark qr($A4)
# BenchmarkTools.Trial: 414 samples with 1 evaluation.
#  Range (min … max):   9.350 ms … 17.288 ms  ┊ GC (min … max):  0.00% … 37.45%
#  Time  (median):     12.411 ms              ┊ GC (median):     9.47%
#  Time  (mean ± σ):   12.068 ms ±  1.869 ms  ┊ GC (mean ± σ):  11.46% ± 10.23%

#        ▃█▄▇▄                  ▁ ▄▂▁▃▃                          
#   ▃▃▄▆▇█████▆▄▄▂▂▃▃▅▃▂▂▃▁▂▃▃▁▆█▅█████▇▆▅▁▄▂▃▃▁▁▃▃▃▃▄▂▃▄▃▃▂▃▃▂ ▄
#   9.35 ms         Histogram: frequency by time        16.5 ms <

#  Memory estimate: 25.44 MiB, allocs estimate: 167.

b = randn(n)

@benchmark $A1 \ $b
# BenchmarkTools.Trial: 5771 samples with 1 evaluation.
#  Range (min … max):  820.926 μs …   4.211 ms  ┊ GC (min … max): 0.00% … 75.87%
#  Time  (median):     836.495 μs               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   864.003 μs ± 253.387 μs  ┊ GC (mean ± σ):  2.46% ±  6.52%

#   ▅█▇▅▄▂                                                        ▁
#   ███████▇▅▆▆▆▇▇▅▆▄▅▁▁▁▃▃▁▃▁▁▁▁▁▃▁▁▁▁▃▁▁▁▄▃▃▃▁▃▃▃▃▁▁▁▃▁▁▁▁▃▄▄▄▄ █
#   821 μs        Histogram: log(frequency) by time       1.31 ms <

#  Memory estimate: 344.67 KiB, allocs estimate: 16.

@benchmark $A2 \ $b
# BenchmarkTools.Trial: 685 samples with 1 evaluation.
#  Range (min … max):  6.424 ms … 12.427 ms  ┊ GC (min … max): 0.00% … 21.60%
#  Time  (median):     6.797 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   7.278 ms ±  1.040 ms  ┊ GC (mean ± σ):  5.60% ±  9.66%

#      ▃▇█▁                                                     
#   ▃▄▆████▆▄▃▃▃▃▃▃▂▃▃▂▂▂▁▂▂▂▁▂▁▁▁▂▁▁▂▁▁▁▁▂▃▃▃▄▃▄▃▃▂▃▂▂▂▁▂▁▂▂▂ ▃
#   6.42 ms        Histogram: frequency by time        10.2 ms <

#  Memory estimate: 7.64 MiB, allocs estimate: 4.

@benchmark $A3 \ $b
# BenchmarkTools.Trial: 742 samples with 1 evaluation.
#  Range (min … max):  3.821 ms … 19.027 ms  ┊ GC (min … max):  0.00% … 54.32%
#  Time  (median):     5.582 ms              ┊ GC (median):     0.00%
#  Time  (mean ± σ):   6.733 ms ±  3.499 ms  ┊ GC (mean ± σ):  20.65% ± 22.34%

#   █▇▁▂ ▁▁▁ ▁ ▁ ▁▂▆▅▂                          ▁               
#   ████▇█████████████▄▄▅▄▄▁▄▁▁▁▁▁▄▁▁▁▁▁▁▁▇▇██▆████▇██▆█▇▆▆█▁▅ █
#   3.82 ms      Histogram: log(frequency) by time     15.9 ms <

#  Memory estimate: 7.55 MiB, allocs estimate: 78253.

@benchmark $A4 \ $b
# BenchmarkTools.Trial: 2510 samples with 1 evaluation.
#  Range (min … max):  1.753 ms …   5.893 ms  ┊ GC (min … max): 0.00% … 65.44%
#  Time  (median):     1.843 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   1.988 ms ± 624.193 μs  ┊ GC (mean ± σ):  6.63% ± 12.28%

#   ▆█▇▄▁                                                       ▁
#   ██████▆▃▅▄▁▁▁▅▁▃▁▁▁▁▁▁▁▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▇███▇▇▆ █
#   1.75 ms      Histogram: log(frequency) by time      4.78 ms <

#  Memory estimate: 2.21 MiB, allocs estimate: 82.
```

## Public API

For documentation of the public API, please use the REPL help mode.

```
AlmostBandedMatrix
bandpart
fillpart
exclusive_bandpart
finish_part_setindex!
```

## Some Considerations

1. Not all operations are meant to be fast. Currently only the following operations are
   fast or planned to be fast:

    + [x] Matrix Vector Multiply
    + [x] Matrix Matrix Multiply
    + [x] QR Factorization
    + [x] UpperTriangular ldiv!
    + [ ] SubArrays

2. Broadcasting won't be fast and will materialize the array into a standard Julia Array.

3. ArrayInferface.jl API

    + [x] `fast_scalar_indexing`
    + [x] `qr_instance`

4. No support for `fast_scalar_indexing == false` arrays.

5. `LinearSolve.jl` has proper dispatches setup.
