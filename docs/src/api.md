# [API](@id API)

The following names make up the public API of FastAlmostBandedMatrices.jl.

## Matrix Type

```@docs
AlmostBandedMatrix
```

## Banded Construction Helpers

`FastAlmostBandedMatrices` reexports the `BandedMatrices` constructor helpers needed by the
documented construction examples: `BandedMatrix` for explicit banded storage and `brand`
for random banded test matrices.

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
