using SafeTestsets

@safetestset "Reexported BandedMatrices API" begin
    using FastAlmostBandedMatrices
    import BandedMatrices

    # Kept in sync with the reexport `export` block in src/FastAlmostBandedMatrices.jl,
    # `REEXPORTED_API` in test/qa/qa.jl, and the reexport section of docs/src/api.md.
    reexports = (
        :Band, :BandError, :BandRange, :BandedMatrix, :band, :bandrange, :bandwidth,
        :bandwidths, :brand, :brandn, :colrange, :rowrange,
    )

    exported = Set(names(FastAlmostBandedMatrices))
    for name in reexports
        # Exported, in scope from a bare `using FastAlmostBandedMatrices`, and still the
        # upstream binding rather than a copy.
        @test name in exported
        @test isdefined(@__MODULE__, name)
        @test getfield(FastAlmostBandedMatrices, name) === getfield(BandedMatrices, name)
    end

    # Whole dependencies are not reexported: BandedMatrices' own reexports (FillArrays)
    # and unrelated names stay out.
    for name in (:Fill, :Ones, :Zeros, :Eye, :symrcm, :BandedMatrices)
        @test name ∉ exported
    end

    # The same list lives in three places: the `export` block in src (checked above via
    # `names`), `REEXPORTED_API` in test/qa/qa.jl, and the "Reexported from
    # BandedMatrices.jl" section of docs/src/api.md. They must not drift apart.
    qa_file = joinpath(@__DIR__, "qa", "qa.jl")
    if isfile(qa_file)
        block = match(r"REEXPORTED_API = \((.*?)\)"s, read(qa_file, String))
        @test block !== nothing
        declared = Set(Symbol(m[1]) for m in eachmatch(r":(\w+)", block[1]))
        @test declared == Set(reexports)
    end

    docs_file = joinpath(@__DIR__, "..", "docs", "src", "api.md")
    if isfile(docs_file)
        section = match(
            r"## Reexported from BandedMatrices\.jl(.*?)\n\nThese names are owned"s,
            read(docs_file, String),
        )
        @test section !== nothing
        bullets = join(filter(startswith("  - "), split(section[1], '\n')), '\n')
        documented = Set(Symbol(m[1]) for m in eachmatch(r"`(\w+)`", bullets))
        @test documented == Set(reexports)
    end

    # The documented construction path works with `using FastAlmostBandedMatrices` alone.
    A = AlmostBandedMatrix(brand(Float64, 10, 10, 3, 2), rand(Float64, 2, 10))
    B = AlmostBandedMatrix(BandedMatrix(fill(1.0, 10, 10), (3, 2)), rand(Float64, 2, 10))

    @test A isa AlmostBandedMatrix
    @test B isa AlmostBandedMatrix
    @test bandpart(A) isa BandedMatrix
    @test bandwidths(bandpart(A)) == (3, 2)
    @test almostbandwidths(A) == (3, 2)
    @test almostbandedrank(A) == 2
end

@safetestset "Constructors" begin
    using FastAlmostBandedMatrices

    A = AlmostBandedMatrix{Float64}(undef, (10, 11), (2, 1), 2)
    A[1, 1] = 2
    @test A[1, 1] == 2.0
    A[4, 1] = 0
    @test A[4, 1] == 0.0
    @test_throws BandError A[4, 1] = 2
    A[1, 3] = 5
    @test A[1, 3] == 5.0

    @test almostbandwidths(A) == (2, 1)
    @test almostbandedrank(A) == 2
end

@safetestset "similar" begin
    using FastAlmostBandedMatrices

    A = AlmostBandedMatrix(brand(Float64, 10, 10, 2, 1), rand(Float64, 2, 10))

    @test similar(A) isa AlmostBandedMatrix
    @test similar(A, Float32) isa AlmostBandedMatrix{Float32}

    fallback = similar(A, Float32, 10, 10)
    @test fallback isa Matrix{Float32}
    @test size(fallback) == size(A)
end

@safetestset "Copy" begin
    using FastAlmostBandedMatrices

    n = 5
    m = 2

    A1 = AlmostBandedMatrix(brand(Float64, n, n, m + 1, m), rand(Float64, m, n))
    A2 = copy(A1)

    @test !(A2 isa Matrix)
    @test A1 == A2

    A2 = deepcopy(A1)

    @test !(A2 isa Matrix)
    @test A1 == A2
end

@safetestset "QR" begin
    using LinearAlgebra, FastAlmostBandedMatrices
    import MatrixFactorizations: QRPackedQ

    n = 80
    A = AlmostBandedMatrix(BandedMatrix(fill(2.0, n, n), (1, 1)), fill(3.0, 2, n))
    A[band(0)] .+= 1:n
    Ã = deepcopy(A)
    B, L = bandpart(A), fillpart(A)

    F = qr(A)
    @test F.Q isa LinearAlgebra.QRPackedQ{Float64, <:BandedMatrix}
    @test F.R isa UpperTriangular{Float64, <:SubArray{Float64, 2, <:AlmostBandedMatrix}}
    @test F.Q' * A ≈ F.R
    @test A == Ã

    @inferred qr(A)

    b = randn(n)
    @test A \ b ≈ Matrix(A) \ b
    @test all(A \ b .=== F \ b)
    @test all(A \ b .=== F.R \ (F.Q' * b))
    Q̃ = QRPackedQ(F.factors, F.τ)
    @test Matrix(Q̃) ≈ Matrix(F.Q)
    @test lmul!(Q̃, copy(b)) ≈ lmul!(F.Q, copy(b)) ≈ Matrix(F.Q) * b
    @test lmul!(Q̃', copy(b)) ≈ lmul!(F.Q', copy(b)) ≈ Matrix(F.Q)' * b
end

@safetestset "Triangular" begin
    using LinearAlgebra, ArrayLayouts, FastAlmostBandedMatrices
    import FastAlmostBandedMatrices: AlmostBandedLayout

    n = 80
    A = AlmostBandedMatrix(BandedMatrix(fill(2.0, n, n), (1, 1)), fill(3.0, 1, n))
    b = randn(n)
    @test MemoryLayout(UpperTriangular(A)) ==
        TriangularLayout{'U', 'N', AlmostBandedLayout}()
    @test_broken UpperTriangular(Matrix(A)) \ b ≈ UpperTriangular(A) \ b
    @test_broken UnitUpperTriangular(Matrix(A)) \ b ≈ UnitUpperTriangular(A) \ b
    @test LowerTriangular(Matrix(A)) \ b ≈ LowerTriangular(A) \ b
    @test UnitLowerTriangular(Matrix(A)) \ b ≈ UnitLowerTriangular(A) \ b
end

# https://github.com/SciML/FastAlmostBandedMatrices.jl/issues/19
@safetestset "fill! on sparse array with BigFloat" begin
    using FastAlmostBandedMatrices, SparseArrays

    A = sparse([1, 2], [1, 5], big.([1.0, 1.0]))
    A1 = AlmostBandedMatrix(brand(BigFloat, 5, 5, 1, 1), A)
    fill!(A1, BigFloat(0.0))
    @test length(A1.fill.nzval) == 2
end
