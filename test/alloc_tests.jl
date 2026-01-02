using AllocCheck
using BenchmarkTools
using FastAlmostBandedMatrices
using FastAlmostBandedMatrices: DisjointRange
using ArrayLayouts: colsupport, rowsupport
using Test

@testset "Allocation Tests" begin
    @testset "DisjointRange - Zero Allocations" begin
        # Test that DisjointRange operations don't allocate
        r1 = Base.OneTo(5)
        r2 = 10:15
        dr = DisjointRange(r1, r2)

        # Test length
        allocs = @allocated length(dr)
        @test allocs == 0

        # Test getindex
        allocs = @allocated dr[3]
        @test allocs == 0

        allocs = @allocated dr[8]
        @test allocs == 0

        # Test first/last
        allocs = @allocated first(dr)
        @test allocs == 0

        allocs = @allocated last(dr)
        @test allocs == 0

        # Test iteration (after warmup)
        sum_test = 0
        for x in dr
            sum_test += x
        end
        allocs = @allocated begin
            s = 0
            for x in dr
                s += x
            end
            s
        end
        @test allocs == 0
    end

    @testset "colsupport - Zero Allocations" begin
        n = 100
        m = 2
        B = brand(Float64, n, n, m + 1, m)
        F = rand(Float64, m, n)
        A = AlmostBandedMatrix(B, F)

        # Warmup
        colsupport(A, 5)
        colsupport(A, 50)

        # Test colsupport for j <= l+u (should return OneTo, no allocation)
        allocs = @allocated colsupport(A, 5)
        @test allocs == 0

        # Test colsupport for j > l+u (now returns DisjointRange instead of vcat)
        allocs = @allocated colsupport(A, 50)
        @test allocs == 0
    end

    @testset "rowsupport - Zero Allocations" begin
        n = 100
        m = 2
        B = brand(Float64, n, n, m + 1, m)
        F = rand(Float64, m, n)
        A = AlmostBandedMatrix(B, F)

        # Warmup
        rowsupport(A, 1)
        rowsupport(A, 50)

        # Test rowsupport (always returns UnitRange, no allocation)
        allocs = @allocated rowsupport(A, 1)
        @test allocs == 0

        allocs = @allocated rowsupport(A, 50)
        @test allocs == 0
    end

    @testset "getindex/setindex! - Zero Allocations" begin
        n = 100
        m = 2
        B = brand(Float64, n, n, m + 1, m)
        F = rand(Float64, m, n)
        A = AlmostBandedMatrix(B, F)

        # Warmup
        _ = A[50, 50]
        A[50, 50] = 1.0

        # Test getindex
        allocs = @allocated A[50, 50]
        @test allocs == 0

        # Test setindex! in band part
        allocs = @allocated A[50, 50] = 2.0
        @test allocs == 0

        # Test setindex! in fill part
        allocs = @allocated A[1, 50] = 3.0
        @test allocs == 0

        # Test setindex! in overlapping part
        allocs = @allocated A[1, 1] = 4.0
        @test allocs == 0
    end

    @testset "bandpart/fillpart - Zero Allocations" begin
        n = 100
        m = 2
        B = brand(Float64, n, n, m + 1, m)
        F = rand(Float64, m, n)
        A = AlmostBandedMatrix(B, F)

        # Warmup
        bandpart(A)
        fillpart(A)

        # Test bandpart
        allocs = @allocated bandpart(A)
        @test allocs == 0

        # Test fillpart
        allocs = @allocated fillpart(A)
        @test allocs == 0
    end
end
