using SafeTestsets, Test

@testset "FastAlmostBandedMatrices" begin
    @safetestset "Constructors" begin
        using FastAlmostBandedMatrices

        A = AlmostBandedMatrix{Float64}(undef, (10, 11), (2, 1), 2)
        A[1, 1] = 2
        @test A[1, 1] == 2.0
        A[4, 1] = 0
        @test A[4, 1] == 0.0
        @test_throws BandError A[4, 1]=2
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
        Ã = deepcopy(A)
        B, L = bandpart(A), fillpart(A)

        F = qr(A)
        @test F.Q isa LinearAlgebra.QRPackedQ{Float64, <:BandedMatrix}
        @test F.R isa UpperTriangular{Float64, <:SubArray{Float64, 2, <:AlmostBandedMatrix}}
        @test F.Q' * A ≈ F.R
        @test A == Ã

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

    @safetestset "Aqua Q/A" begin
        using Aqua, FastAlmostBandedMatrices

        Aqua.test_all(FastAlmostBandedMatrices; ambiguities = false)
    end
end
