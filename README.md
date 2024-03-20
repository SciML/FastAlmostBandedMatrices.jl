# FastAlmostBandedMatrices.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)

[![CI](https://github.com/SciML/FastAlmostBandedMatrices.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/SciML/FastAlmostBandedMatrices.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/SciML/FastAlmostBandedMatrices.jl/branch/main/graph/badge.svg?)](https://codecov.io/gh/SciML/FastAlmostBandedMatrices.jl)
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

### QR Factorization & Linear Solve

<details>
 <summary>Click me!</summary>
<p>

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
# BenchmarkTools.Trial: 6959 samples with 1 evaluation.
#  Range (min … max):  680.420 μs …   2.600 ms  ┊ GC (min … max): 0.00% … 70.11%
#  Time  (median):     707.540 μs               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   716.433 μs ± 100.758 μs  ┊ GC (mean ± σ):  1.05% ±  4.81%

#            ▃▃▅▅██▃▂                                              
#   ▂▂▂▃▃▃▂▃▇█████████▆▅▄▄▃▃▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▁▁▁▂▁▁▂▂▂▂▂▂▁▂▂▂▂▂▁▂ ▃
#   680 μs           Histogram: frequency by time          802 μs <

#  Memory estimate: 320.70 KiB, allocs estimate: 11.

@benchmark qr($A2)
# BenchmarkTools.Trial: 173 samples with 1 evaluation.
#  Range (min … max):  23.543 ms … 39.269 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     29.162 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   28.962 ms ±  2.984 ms  ┊ GC (mean ± σ):  0.75% ± 2.29%

#            ▂▁▃  █       ▂▄▄▂▂ ▁                                
#   ▅▅▁▃▃▃▅▁▇███▇▆██▅▄▃▇▄▃███████▅▄▅▄▃▃▁▅▃▁▃▃▄▅▁▃▃▄▃▁▁▃▁▁▁▁▁▁▁▃ ▃
#   23.5 ms         Histogram: frequency by time        38.4 ms <

#  Memory estimate: 8.18 MiB, allocs estimate: 6.

@benchmark qr($A3)
# BenchmarkTools.Trial: 1452 samples with 1 evaluation.
#  Range (min … max):  2.975 ms …   7.381 ms  ┊ GC (min … max): 0.00% … 21.67%
#  Time  (median):     3.206 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   3.440 ms ± 625.169 μs  ┊ GC (mean ± σ):  6.43% ± 11.76%

#     ▂▆██▇▅▂                                ▁▁  ▂▁▁▁            
#   ██████████▆▅▅▄▄▄▄▁▁▄▁▄▄▁▁▁▁▄▁▄▄▄▁▄▁▁▁▁▁▆▇████████▇▅▇▅▅▅▁▄▁▅ █
#   2.97 ms      Histogram: log(frequency) by time      5.52 ms <

#  Memory estimate: 6.32 MiB, allocs estimate: 65108.

@benchmark qr($A4)
# BenchmarkTools.Trial: 362 samples with 1 evaluation.
#  Range (min … max):   8.413 ms … 104.416 ms  ┊ GC (min … max): 0.00% … 1.04%
#  Time  (median):     10.826 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   13.788 ms ±   9.973 ms  ┊ GC (mean ± σ):  1.20% ± 3.47%

#    ▄█                                                           
#   ▄███▆▄▆▅▅▅▄▃▂▃▃▂▂▂▃▃▃▄▃▁▂▁▂▁▂▁▃█▇█▅▆▆▅▃▄▃▃▂▃▂▃▁▃▂▃▃▂▃▂▂▂▃▃▁▃ ▃
#   8.41 ms         Histogram: frequency by time         22.3 ms <

#  Memory estimate: 25.44 MiB, allocs estimate: 169.

b = randn(n)

@benchmark $A1 \ $b
# BenchmarkTools.Trial: 5945 samples with 1 evaluation.
#  Range (min … max):  797.407 μs …   2.972 ms  ┊ GC (min … max): 0.00% … 68.87%
#  Time  (median):     828.066 μs               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   838.834 μs ± 124.996 μs  ┊ GC (mean ± σ):  1.27% ±  5.33%

#                 ▁▂▂▇█▇▃▂▁                                        
#   ▂▂▃▄▅▆▇█▇▆█▇███████████▇█▇▆▆▅▄▄▃▃▃▃▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂▂▁▁▂▁▂▂▁▁▂▂ ▄
#   797 μs           Histogram: frequency by time          901 μs <

#  Memory estimate: 367.95 KiB, allocs estimate: 431.

@benchmark $A2 \ $b
# BenchmarkTools.Trial: 586 samples with 1 evaluation.
#  Range (min … max):  7.682 ms …  15.557 ms  ┊ GC (min … max): 0.00% … 40.03%
#  Time  (median):     8.305 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   8.515 ms ± 679.543 μs  ┊ GC (mean ± σ):  1.86% ±  3.97%

#        ▂█▃▂▅▆▄     ▁                                           
#   ▄▆▆▅▆█████████▅▇▆█▄▄▃▄▆▆▅▆▄▃▅▃▄▃▄▂▃▁▃▃▃▃▁▂▂▃▁▃▃▁▂▁▁▁▃▁▁▁▁▃▂ ▃
#   7.68 ms         Histogram: frequency by time          11 ms <

#  Memory estimate: 7.64 MiB, allocs estimate: 4.

@benchmark $A3 \ $b
# BenchmarkTools.Trial: 1212 samples with 1 evaluation.
#  Range (min … max):  3.458 ms …   7.806 ms  ┊ GC (min … max): 0.00% … 30.80%
#  Time  (median):     3.815 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   4.118 ms ± 707.897 μs  ┊ GC (mean ± σ):  7.61% ± 12.33%

#       ▂█▇▆▆                                                    
#   ▂▃▄▅██████▆▄▃▃▃▂▂▂▂▂▂▂▂▁▂▁▂▁▂▂▁▁▁▂▂▂▂▃▃▄▄▄▃▃▃▃▃▃▃▂▂▂▃▂▁▂▂▂▂ ▃
#   3.46 ms         Histogram: frequency by time        6.17 ms <

#  Memory estimate: 7.55 MiB, allocs estimate: 78253.

@benchmark $A4 \ $b
# BenchmarkTools.Trial: 2879 samples with 1 evaluation.
#  Range (min … max):  1.627 ms …   3.120 ms  ┊ GC (min … max): 0.00% … 36.35%
#  Time  (median):     1.675 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   1.734 ms ± 243.465 μs  ┊ GC (mean ± σ):  2.76% ±  7.52%

#   ▂▇█▆▃                                                   ▁▁  ▁
#   ██████▇█▆▄▅▃▄▁▁▁▁▃▁▁▁▁▁▁▃▄▁▁▃▁▃▁▁▁▃▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▇███ █
#   1.63 ms      Histogram: log(frequency) by time      2.81 ms <

#  Memory estimate: 2.21 MiB, allocs estimate: 82.
```

</p>
</details>

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
    + [x] SubArrays

2. Broadcasting won't be fast and will materialize the array into a standard Julia Array.

3. ArrayInferface.jl API

    + [x] `fast_scalar_indexing`
    + [x] `qr_instance`

4. No support for `fast_scalar_indexing == false` arrays.

5. `LinearSolve.jl` has proper dispatches setup.
