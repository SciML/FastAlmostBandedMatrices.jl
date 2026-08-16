# [API](@id API)

The following names make up the public API of FastAlmostBandedMatrices.jl.

## Matrix Type

```@docs
AlmostBandedMatrix
```

## Banded Construction Helpers

`BandedMatrix` and `brand` are reexported because they are the public constructors used by
the documented `AlmostBandedMatrix` construction workflow. Other banded-matrix utilities
should be imported from `BandedMatrices` directly.

```@docs
BandedMatrix
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
