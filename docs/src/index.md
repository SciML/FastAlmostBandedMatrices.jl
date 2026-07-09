# FastAlmostBandedMatrices.jl

A fast implementation of Almost Banded Matrices in Julia. Primarily developed for use in
[BoundaryValueDiffEq.jl](https://github.com/SciML/BoundaryValueDiffEq.jl).

An almost banded matrix is a banded matrix together with a low-rank "fill" part occupying its
leading rows. Storing and factorizing the two parts separately allows QR factorizations and
linear solves that are much faster than treating the matrix as dense or sparse.

[SemiseparableMatrices.jl](https://github.com/JuliaLinearAlgebra/SemiseparableMatrices.jl)
provides an alternate implementation of almost banded matrices, but that had considerable
overhead for this use-case. This implementation borrows a lot of details from that
repository.

## Installation

```julia
using Pkg
Pkg.add("FastAlmostBandedMatrices")
```

## Basic Construction and Usage

```julia
using FastAlmostBandedMatrices, LinearAlgebra

m = 2  # Fill rank
n = 10 # Matrix dimension

# Create an almost banded matrix.
# brand(T, m, n, lower_bandwidth, upper_bandwidth) creates a random banded matrix.
A = AlmostBandedMatrix(brand(Float64, n, n, m + 1, m), rand(Float64, m, n))

# Query matrix properties
almostbandwidths(A)  # Returns (3, 2) - the bandwidths of the banded part
almostbandedrank(A)  # Returns 2 - the rank of the fill part

# Access the parts
B = bandpart(A)       # Get the banded part
F = fillpart(A)       # Get the fill part (low-rank correction)

# Solve linear systems efficiently
b = rand(n)
x = A \ b             # Uses a specialized QR factorization

# QR factorization
fact = qr(A)
Q, R = fact.Q, fact.R
```

## Public API

See the [API](@ref) page for the documented public interface.
