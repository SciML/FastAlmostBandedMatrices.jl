using FastAlmostBandedMatrices, BenchmarkTools
using StableRNGs, LinearAlgebra

const SUITE = BenchmarkGroup()
const rng = StableRNG(123)

# AlmostBandedMatrix = banded part + full "fill" rows at the bottom
n = 300
m = 4
banded_part = BandedMatrix(rand(rng, n, n), (3, 2))
fill_part = rand(rng, m, n)
A = AlmostBandedMatrix(banded_part, fill_part)
x = rand(rng, n)
y = zeros(n)

A_dense = Matrix(A)

# =============================================================================
# Construction
# =============================================================================

SUITE["construct"] = BenchmarkGroup()

SUITE["construct"]["almost_banded"] = @benchmarkable AlmostBandedMatrix(
    $banded_part, $fill_part
)
SUITE["construct"]["banded"] = @benchmarkable BandedMatrix(
    $(rand(rng, n, n)), (3, 2)
)

# =============================================================================
# Operations
# =============================================================================

SUITE["ops"] = BenchmarkGroup()

SUITE["ops"]["matvec"] = @benchmarkable $A * $x
SUITE["ops"]["mul!"] = @benchmarkable mul!($y, $A, $x)
SUITE["ops"]["to_dense"] = @benchmarkable Matrix($A)
SUITE["ops"]["lu"] = @benchmarkable lu($A)
SUITE["ops"]["ldiv"] = @benchmarkable $A_dense \ $x
