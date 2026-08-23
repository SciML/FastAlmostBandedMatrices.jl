# [API](@id API)

The following names make up the public API of FastAlmostBandedMatrices.jl.

## Matrix Type

```@docs
AlmostBandedMatrix
```

## Accessing the Parts

```@docs
bandpart
fillpart
exclusive_bandpart
```

## Querying Structure

```@docs
almostbandwidths
almostbandedrank
```

## Mutation Helpers

```@docs
finish_part_setindex!
```

## Reexported from BandedMatrices.jl

The documented constructor is `AlmostBandedMatrix(bands::BandedMatrix, fill)`, and the
examples on the home page build `bands` with `brand`. So `using FastAlmostBandedMatrices`
also brings the
[BandedMatrices.jl](https://github.com/JuliaLinearAlgebra/BandedMatrices.jl) names needed to
construct, populate and query that banded argument:

| Purpose | Names |
|:--------|:------|
| Construction | `BandedMatrix`, `brand`, `brandn` |
| Band indexing | `Band`, `BandRange`, `band`, `bandrange` |
| Structure queries | `bandwidth`, `bandwidths`, `colrange`, `rowrange` |
| Errors | `BandError` |

These names are owned and documented by BandedMatrices.jl; see its documentation for their
full behaviour. Nothing else from BandedMatrices.jl is reexported, and neither are the
names it reexports from its own dependencies (for example the FillArrays.jl names `Fill`,
`Ones`, `Zeros` and `Eye`). Use `using BandedMatrices` or `using FillArrays` for those.
