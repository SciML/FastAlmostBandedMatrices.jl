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

`AlmostBandedMatrix(bands::BandedMatrix, fill)` is the documented constructor, and the
examples on the home page build `bands` with `brand`. So `using FastAlmostBandedMatrices`
also brings in the [BandedMatrices.jl](https://juliaLinearAlgebra.github.io/BandedMatrices.jl/stable/)
names needed to build that banded argument, populate it and query it, so they do not have
to be imported separately:

  - Building the bands: `BandedMatrix`, `brand`, `brandn`
  - Populating them: `Band`, `BandRange`, `band`, `bandrange`
  - Querying them: `bandwidth`, `bandwidths`, `colrange`, `rowrange`
  - Errors: `BandError`

These names are defined and documented by
[BandedMatrices.jl](https://juliaLinearAlgebra.github.io/BandedMatrices.jl/stable/).

Anything else from BandedMatrices.jl must be imported from BandedMatrices.jl directly.
In particular, the following are deliberately **not** reexported:

  - The [FillArrays.jl](https://juliaarrays.github.io/FillArrays.jl/stable/) names
    `Fill`, `Ones`, `Zeros` and `Eye`. BandedMatrices.jl reexports those from a
    dependency of its own, so passing them along here would be reexporting a whole
    dependency chain; none of them is needed to build an `AlmostBandedMatrix`. Use
    `using FillArrays`.
  - The `BandedMatrices` module binding itself. Use `import BandedMatrices` if you want
    the qualified form.
  - `symrcm`, a sparse-matrix reordering unrelated to this package.

The `Band` and `BandError` types are defined by
[BandedMatrices.jl](https://juliaLinearAlgebra.github.io/BandedMatrices.jl/stable/), and
FastAlmostBandedMatrices.jl documents their public bindings here because they are part of
its documented construction workflow. Every other reexported name above carries its
upstream docstring, reachable from the REPL help mode.

```@docs
Band
BandError
```

The list above is kept in sync with the reexport `export` block in
`src/FastAlmostBandedMatrices.jl` and with `REEXPORTED_API` in `test/qa/qa.jl`, and the
Core test suite asserts that all three agree.
