using AllocCheck
using FastAlmostBandedMatrices
using FastAlmostBandedMatrices: DisjointRange
using ArrayLayouts: colsupport, rowsupport
using Test

# Allocation freedom is verified with `AllocCheck.check_allocs`, which statically
# proves the absence of any allocating code path for the given argument types.
# This is immune to the first-call compilation and global-scope boxing noise that
# `@allocated` picks up on newer Julia versions, so it directly tests the intended
# invariant (the operation never allocates) rather than a single runtime sample.

getidx(A, i...) = A[i...]
setidx!(A, v, i...) = (A[i...] = v; nothing)
function sumiter(d)
    s = zero(eltype(d))
    for x in d
        s += x
    end
    return s
end

@testset "Allocation Tests" begin
    @testset "DisjointRange - Zero Allocations" begin
        r1 = Base.OneTo(5)
        r2 = 10:15
        dr = DisjointRange(r1, r2)

        @test isempty(check_allocs(length, (typeof(dr),)))
        @test isempty(check_allocs(getidx, (typeof(dr), Int)))
        @test isempty(check_allocs(first, (typeof(dr),)))
        @test isempty(check_allocs(last, (typeof(dr),)))
        @test isempty(check_allocs(sumiter, (typeof(dr),)))
    end

    @testset "colsupport - Zero Allocations" begin
        n = 100
        m = 2
        B = brand(Float64, n, n, m + 1, m)
        F = rand(Float64, m, n)
        A = AlmostBandedMatrix(B, F)

        # colsupport returns OneTo for j <= l+u and a DisjointRange otherwise; both
        # branches must be allocation-free.
        @test isempty(check_allocs(colsupport, (typeof(A), Int)))
    end

    @testset "rowsupport - Zero Allocations" begin
        n = 100
        m = 2
        B = brand(Float64, n, n, m + 1, m)
        F = rand(Float64, m, n)
        A = AlmostBandedMatrix(B, F)

        @test isempty(check_allocs(rowsupport, (typeof(A), Int)))
    end

    @testset "getindex/setindex! - Zero Allocations" begin
        n = 100
        m = 2
        B = brand(Float64, n, n, m + 1, m)
        F = rand(Float64, m, n)
        A = AlmostBandedMatrix(B, F)

        @test isempty(check_allocs(getidx, (typeof(A), Int, Int)))
        @test isempty(check_allocs(setidx!, (typeof(A), Float64, Int, Int)))
    end

    @testset "bandpart/fillpart - Zero Allocations" begin
        n = 100
        m = 2
        B = brand(Float64, n, n, m + 1, m)
        F = rand(Float64, m, n)
        A = AlmostBandedMatrix(B, F)

        @test isempty(check_allocs(bandpart, (typeof(A),)))
        @test isempty(check_allocs(fillpart, (typeof(A),)))
    end
end
