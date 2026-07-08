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

# On aarch64 (e.g. Apple Silicon) the per-task GC stack pointer is fetched through
# `jl_get_pgcstack_fallback` rather than a TLS/register fast path. AllocCheck already
# treats the fast-path variants `jl_get_pgcstack`/`jl_get_pgcstack_static` as
# non-allocating, but its whitelist omits the `_fallback` variant, so it conservatively
# reports it as an `AllocatingRuntimeCall`. That call does not allocate, so we drop it
# before asserting; every other (genuinely allocating) result is still caught.
real_allocs(results) = filter(results) do r
    !(r isa AllocCheck.AllocatingRuntimeCall && r.name == "jl_get_pgcstack_fallback")
end

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

        @test isempty(real_allocs(check_allocs(length, (typeof(dr),))))
        @test isempty(real_allocs(check_allocs(getidx, (typeof(dr), Int))))
        @test isempty(real_allocs(check_allocs(first, (typeof(dr),))))
        @test isempty(real_allocs(check_allocs(last, (typeof(dr),))))
        @test isempty(real_allocs(check_allocs(sumiter, (typeof(dr),))))
    end

    @testset "colsupport - Zero Allocations" begin
        n = 100
        m = 2
        B = brand(Float64, n, n, m + 1, m)
        F = rand(Float64, m, n)
        A = AlmostBandedMatrix(B, F)

        # colsupport returns OneTo for j <= l+u and a DisjointRange otherwise; both
        # branches must be allocation-free.
        @test isempty(real_allocs(check_allocs(colsupport, (typeof(A), Int))))
    end

    @testset "rowsupport - Zero Allocations" begin
        n = 100
        m = 2
        B = brand(Float64, n, n, m + 1, m)
        F = rand(Float64, m, n)
        A = AlmostBandedMatrix(B, F)

        @test isempty(real_allocs(check_allocs(rowsupport, (typeof(A), Int))))
    end

    @testset "getindex/setindex! - Zero Allocations" begin
        n = 100
        m = 2
        B = brand(Float64, n, n, m + 1, m)
        F = rand(Float64, m, n)
        A = AlmostBandedMatrix(B, F)

        @test isempty(real_allocs(check_allocs(getidx, (typeof(A), Int, Int))))
        @test isempty(real_allocs(check_allocs(setidx!, (typeof(A), Float64, Int, Int))))
    end

    @testset "bandpart/fillpart - Zero Allocations" begin
        n = 100
        m = 2
        B = brand(Float64, n, n, m + 1, m)
        F = rand(Float64, m, n)
        A = AlmostBandedMatrix(B, F)

        @test isempty(real_allocs(check_allocs(bandpart, (typeof(A),))))
        @test isempty(real_allocs(check_allocs(fillpart, (typeof(A),))))
    end
end
