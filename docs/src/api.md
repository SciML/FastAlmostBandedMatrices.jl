# [API](@id API)

The following names make up the public API of FastAlmostBandedMatrices.jl.

## Matrix Type

```@docs
AlmostBandedMatrix
```

## Banded Construction Helper

`brand` is reexported because it is used by the documented `AlmostBandedMatrix` construction
workflow. Other banded-matrix types and utilities should be imported from `BandedMatrices`
directly.

```@docs
brand
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
