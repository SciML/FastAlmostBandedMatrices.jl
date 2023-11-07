# FastAlmostBandedMatrices.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)

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
```

## Benchmarks

```julia
using BenchmarkTools, FastAlmostBandedMatrices, SparseArrays, FillArrays, LinearAlgebra
import SemiseparableMatrices

m = 5
n = 1000

II = zeros(n, m)
II[diagind(II)] .= 1

A1 = AlmostBandedMatrix(rand(Float64, m, n), brand(Float64, n, n, m + 1, m))
A2 = Matrix(A1)
A3 = SemiseparableMatrices.AlmostBandedMatrix(copy(A1.bands),
    SemiseparableMatrices.LowRankMatrix(II, copy(A1.fill)))
A4 = sparse(A2)

@benchmark qr($A1)
# BenchmarkTools.Trial: 2708 samples with 1 evaluation.
#  Range (min … max):  1.015 ms … 8.032 ms  ┊ GC (min … max):  0.00% … 56.96%
#  Time  (median):     1.147 ms             ┊ GC (median):     0.00%
#  Time  (mean ± σ):   1.840 ms ± 1.427 ms  ┊ GC (mean ± σ):  15.26% ± 18.02%

#   ▅█▅▂▁                             ▃▁▂▁           ▂▂        
#   ██████▇▇▇▆▅▆▁▃▁▁▁▁▄▃▁▁▄▁▃▄▃▁▁▁▃▁▃▃█████▇▆▅▆▅▄▃▅▅▃███▆▅▄▆▆ █
#   1.02 ms     Histogram: log(frequency) by time     5.96 ms <

#  Memory estimate: 7.94 MiB, allocs estimate: 16.

@benchmark qr($A2)
# BenchmarkTools.Trial: 172 samples with 1 evaluation.
#  Range (min … max):  24.451 ms … 39.738 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     28.915 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   29.063 ms ±  2.653 ms  ┊ GC (mean ± σ):  0.92% ± 2.65%

#       ▃               █▂  ▂  ▂               ▂                 
#   ▆▄▃▅█▄▃▅▃▄▆▄▅▄▄▄▄▃▅▅██▇▇██▆█▆▅▆▄█▄▇▆▇▆▃▃▃▄▄█▁▃▄▁▃▃▃▄▃▁▁▃▁▃▃ ▃
#   24.5 ms         Histogram: frequency by time          35 ms <

#  Memory estimate: 8.18 MiB, allocs estimate: 6.

@benchmark qr($A3)
# BenchmarkTools.Trial: 1548 samples with 1 evaluation.
#  Range (min … max):  2.791 ms …   6.872 ms  ┊ GC (min … max): 0.00% … 53.51%
#  Time  (median):     2.941 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   3.226 ms ± 791.075 μs  ┊ GC (mean ± σ):  8.59% ± 14.11%

#   ▄▅██▆                                              ▁▁▂▂▁ ▁   
#   █████▇▇▇▆▅▄▆▄▆▄▄▄▄▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▆▇███████▆ █
#   2.79 ms      Histogram: log(frequency) by time      5.53 ms <

#  Memory estimate: 6.32 MiB, allocs estimate: 65108.

@benchmark qr($A4)
# BenchmarkTools.Trial: 441 samples with 1 evaluation.
#  Range (min … max):   8.352 ms … 96.160 ms  ┊ GC (min … max): 0.00% … 1.29%
#  Time  (median):      9.961 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   11.332 ms ±  5.821 ms  ┊ GC (mean ± σ):  3.33% ± 5.91%

#     ▃█▅       ▁                                                
#   ▅████▇▄▃▂▄▇██▇▅▂▂▂▂▂▂▁▃▁▁▁▁▂▁▁▁▁▁▁▁▁▁▁▁▂▃▃▇▆▇▃▄▃▂▂▂▂▁▁▁▁▁▁▂ ▃
#   8.35 ms         Histogram: frequency by time        18.6 ms <

#  Memory estimate: 25.44 MiB, allocs estimate: 169.
```

## Public API

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

    + [ ] Matrix Vector Multiply
    + [ ] Matrix Matrix Multiply
    + [x] QR Factorization
    + [ ] UpperTriangular ldiv!
    + [ ] SubArrays

2. Broadcasting won't be fast and will materialize the array into a standard Julia Array.

3. ArrayInferface.jl API

    + [x] `fast_scalar_indexing`
    + [x] `qr_instance`

4. No support for `fast_scalar_indexing == false` arrays.

5. `LinearSolve.jl` has proper dispatches setup.

6. Internally for QR Factorization, we materialize the entire matrix for fast BLAS calls.
   This would be undesirable for problems where memory is a concern. Use
   `SemiseparableMatrices.jl` for these cases, however, that turns out to be quite slow
   since `1 gemm` call gets expanded into `m gemv` calls where `m` is of the order of the
   dimension of the matrix.
