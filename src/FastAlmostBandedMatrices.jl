module FastAlmostBandedMatrices

import PrecompileTools: @setup_workload, @compile_workload

using ArrayInterface, ArrayLayouts, BandedMatrices, ConcreteStructs, LazyArrays,
      LinearAlgebra, MatrixFactorizations, Reexport

import ArrayLayouts: MemoryLayout, sublayout, sub_materialize, MatLdivVec, materialize!,
                     triangularlayout, triangulardata, zero!, _copyto!, colsupport,
                     rowsupport, _qr, _qr!, _factorize, muladd!
import BandedMatrices: _banded_qr!, bandeddata, banded_qr_lmul!
import LinearAlgebra: ldiv!
import MatrixFactorizations: QR, QRPackedQ, getQ, getR, QRPackedQLayout, AdjQRPackedQLayout

@reexport using BandedMatrices

# ------------------
# DisjointRange - for zero-allocation colsupport
# ------------------

"""
    DisjointRange{T}

A lazy representation of the union of two ranges, supporting iteration and indexing
without heap allocation.
"""
struct DisjointRange{T<:Integer, R1<:AbstractUnitRange{T}, R2<:AbstractUnitRange{T}} <:
       AbstractVector{T}
    r1::R1
    r2::R2
end

Base.size(d::DisjointRange) = (length(d.r1) + length(d.r2),)
Base.length(d::DisjointRange) = length(d.r1) + length(d.r2)

@inline function Base.getindex(d::DisjointRange, i::Integer)
    @boundscheck checkbounds(d, i)
    n1 = length(d.r1)
    if i <= n1
        return @inbounds d.r1[i]
    else
        return @inbounds d.r2[i - n1]
    end
end

Base.IndexStyle(::Type{<:DisjointRange}) = IndexLinear()

@inline function Base.iterate(d::DisjointRange)
    if !isempty(d.r1)
        val, state = iterate(d.r1)
        return val, (1, state)
    elseif !isempty(d.r2)
        val, state = iterate(d.r2)
        return val, (2, state)
    else
        return nothing
    end
end

@inline function Base.iterate(d::DisjointRange, state)
    which, inner_state = state
    if which == 1
        next = iterate(d.r1, inner_state)
        if next !== nothing
            return next[1], (1, next[2])
        else
            # Switch to r2
            if !isempty(d.r2)
                val, new_state = iterate(d.r2)
                return val, (2, new_state)
            else
                return nothing
            end
        end
    else
        next = iterate(d.r2, inner_state)
        if next !== nothing
            return next[1], (2, next[2])
        else
            return nothing
        end
    end
end

Base.first(d::DisjointRange) = isempty(d.r1) ? first(d.r2) : first(d.r1)
Base.last(d::DisjointRange) = isempty(d.r2) ? last(d.r1) : last(d.r2)
Base.minimum(d::DisjointRange) = min(minimum(d.r1), minimum(d.r2))
Base.maximum(d::DisjointRange) = max(maximum(d.r1), maximum(d.r2))

# ------------------
# AlmostBandedMatrix
# ------------------

abstract type AbstractAlmostBandedLayout <: MemoryLayout end
struct AlmostBandedLayout <: AbstractAlmostBandedLayout end

"""
    AlmostBandedMatrix(bands::BandedMatrix, fill)
    AlmostBandedMatrix{T}(bands, fill)
    AlmostBandedMatrix(::UndefInitializer, [::Type{T} = Float64], mn::NTuple{2, Integer},
        lu::NTuple{2, Integer}, rank::Integer)
    AlmostBandedMatrix{T}(::UndefInitializer, mn::NTuple{2, Integer},
        lu::NTuple{2, Integer}, rank::Integer)

An `AlmostBandedMatrix` is a matrix with a `bands` part and a `fill` part. For efficient
operations we store the matrix as a BandedMatrix and another AbstractMatrix with an
overlapping bit.

[3 3 3 2 2 2 2 ... 2 2 2]\\
[3 3 3 3 2 2 2 ... 2 2 2]\\
[0 1 1 1 1 0 0 ... 0 0 0]\\
[0 0 1 1 1 1 0 ... 0 0 0]\\
[.......................]\\
[.......................]\\
[.......................]\\
[0 0 0 0 0 0 0 ... 1 1 0]\\
[0 0 0 0 0 0 0 ... 1 1 1]

where `2`'s are the fill part, and `1`'s are the bands part, and `3`'s are the overlapping
part.
"""
@concrete struct AlmostBandedMatrix{T} <: LayoutMatrix{T}
    bands
    fill
end

function AlmostBandedMatrix(::UndefInitializer, ::Type{T}, mn::NTuple{2, Integer},
        lu::NTuple{2, Integer}, rank::Integer) where {T}
    @assert lu[2] ≥ rank - 1
    @assert rank≥1 "Rank 0 fill array makes it a BandedMatrix."
    bands = BandedMatrix{T}(undef, mn, lu)
    fill = Matrix{T}(undef, rank, mn[2])
    return AlmostBandedMatrix{T}(bands, fill)
end

function AlmostBandedMatrix{T}(::UndefInitializer, mn::NTuple{2, Integer},
        lu::NTuple{2, Integer}, rank::Integer) where {T}
    return AlmostBandedMatrix(undef, T, mn, lu, rank)
end

function AlmostBandedMatrix(::UndefInitializer, mn::NTuple{2, Integer}, lu::NTuple{
        2, Integer}, rank::Integer)
    return AlmostBandedMatrix(undef, Float64, mn, lu, rank)
end

function AlmostBandedMatrix(bands::BandedMatrix, fill::AbstractMatrix)
    @assert size(fill, 2) == size(bands, 2)
    @assert size(fill, 1)≥1 "Rank 0 fill array makes it a BandedMatrix."
    T = promote_type(eltype(fill), eltype(bands))
    @assert bandwidths(bands)[1] ≥ size(fill, 1) - 1
    finish_part_setindex!(bands, fill)
    return AlmostBandedMatrix{T}(bands, fill)
end

MemoryLayout(::Type{<:AlmostBandedMatrix}) = AlmostBandedLayout()

@inline function finish_part_setindex!(A)
    return finish_part_setindex!(bandpart(A), fillpart(A))
end
@inline function finish_part_setindex!(bands, fill)
    # copy `fill` into `bands` in the correct locations
    l, u = bandwidths(bands)
    for i in 1:size(fill, 1), j in max(1, i - l):min(size(bands, 2), i + u)

        @inbounds bands[i, j] = fill[i, j]
    end
    return nothing
end

"""
    bandpart(A::AlmostBandedMatrix)

Banded Part of the AlmostBandedMatrix.
"""
@inline bandpart(A::AlmostBandedMatrix) = A.bands

"""
    fillpart(A::AlmostBandedMatrix)

Fill Part of the AlmostBandedMatrix.
"""
@inline fillpart(A::AlmostBandedMatrix) = A.fill

"""
    exclusive_bandpart(A::AlmostBandedMatrix)

Banded Part of the AlmostBandedMatrix without the overlapping part.
"""
@inline function exclusive_bandpart(A)
    B, F = bandpart(A), fillpart(A)
    return @view(B[(size(F, 1) + 1):end, :])
end

@inline almostbandwidths(_, A) = bandwidths(bandpart(A))
@inline almostbandedrank(_, A) = size(fillpart(A), 1)
@inline almostbandwidths(A) = almostbandwidths(MemoryLayout(typeof(A)), A)
@inline almostbandedrank(A) = almostbandedrank(MemoryLayout(typeof(A)), A)

@inline Base.size(A::AlmostBandedMatrix) = size(A.bands)
@inline Base.eltype(::AlmostBandedMatrix{T}) where {T} = T
@inline Base.IndexStyle(::Type{<:AlmostBandedMatrix}) = IndexCartesian()

# If dims is provided we will construct a Matrix but other invokations should return a
# AlmostBandedMatrix
function Base.similar(A::AlmostBandedMatrix, ::Type{T}) where {T}
    bands = similar(A.bands, T)
    fill = similar(A.fill, T)
    return AlmostBandedMatrix{T}(bands, fill)
end

function Base.fill!(A::AlmostBandedMatrix, v)
    fill!(bandpart(A), v)
    fill!(fillpart(A), v)
    return A
end

@inline function colsupport(::AbstractAlmostBandedLayout, A, j)
    l, u = almostbandwidths(A)
    if j ≤ l + u
        return Base.OneTo(min(maximum(j) + l, size(A, 1)))
    else
        r = almostbandedrank(A)
        sup = colsupport(bandpart(A), j)
        if isempty(sup)
            return Base.OneTo(r)
        else
            # Use DisjointRange to avoid heap allocation from vcat
            return DisjointRange(Base.OneTo(min(r, minimum(sup) - 1)), sup)
        end
    end
end

@inline function rowsupport(::AbstractAlmostBandedLayout, A, k)
    l, _ = almostbandwidths(A)
    if maximum(k) ≤ almostbandedrank(A)
        return max(1, minimum(k) - l):size(A, 2)
    else
        return max(1, minimum(k) - l):min(maximum(k) + l, size(A, 2))
    end
end

@inline function Base.getindex(B::AlmostBandedMatrix, k::Integer, j::Integer)
    if j > k + bandwidth(B.bands, 2) && k ≤ size(B.fill, 1)
        return @inbounds B.fill[k, j]
    else
        return @inbounds B.bands[k, j]
    end
end

function Base.copy(B::AlmostBandedMatrix)
    return AlmostBandedMatrix(copy(bandpart(B)), copy(fillpart(B)))
end

function Base.setindex!(B::AlmostBandedMatrix, v, k::Integer, j::Integer)
    if k ≤ size(B.fill, 1)
        if j > k + bandwidth(B.bands, 2)
            # Only `fill` part
            @inbounds B.fill[k, j] = v
        else
            # Overlapping part
            @inbounds B.bands[k, j] = v
            @inbounds B.fill[k, j] = v
        end
    else
        # Only `band` part
        @inbounds B.bands[k, j] = v
    end
    return
end

function LinearAlgebra.triu!(A::AlmostBandedMatrix)
    triu!(bandpart(A))
    return A
end

# TODO: Support views properly
function sublayout(::AlmostBandedLayout, ::Type{<:Tuple{
        AbstractUnitRange{Int}, AbstractUnitRange{Int}}})
    return AlmostBandedLayout()
end

bandpart(V::SubArray) = view(bandpart(parent(V)), parentindices(V)...)
function fillpart(V::SubArray)
    idx1, idx2 = parentindices(V)
    r = almostbandedrank(parent(V))
    if maximum(idx1) ≤ r
        return view(fillpart(parent(V)), idx1, idx2)
    else
        start = first(idx1)
        if start > r
            return view(fillpart(parent(V)), 0:0, idx2)
        else
            return view(fillpart(parent(V)), start:r, idx2)
        end
    end
    error("Not Implemented!")
end

function almostbandedrank(V::SubArray)
    idx1, _ = parentindices(V)
    r = almostbandedrank(parent(V))
    if maximum(idx1) ≤ r
        return length(idx1)
    else
        return length(first(idx1):r)
    end
end

# ---------------
# Pretty Printing
# ---------------
function _almost_banded_summary(io, B::AlmostBandedMatrix{T}, inds) where {T}
    print(io, Base.dims2string(length.(inds)),
        " AlmostBandedMatrix{$T} with bandwidths $(almostbandwidths(B)) and fill \
          rank $(almostbandedrank(B))")
end
function Base.array_summary(io::IO, B::AlmostBandedMatrix, inds::Tuple{Vararg{Base.OneTo}})
    _almost_banded_summary(io, B, inds)
    print(io, " with data ")
    summary(io, B.bands)
    print(io, " and fill ")
    summary(io, B.fill)
end
function Base.array_summary(io::IO, B::AlmostBandedMatrix, inds)
    _almost_banded_summary(io, B, inds)
    print(io, " with data ")
    summary(io, B.bands)
    print(io, " and fill ")
    summary(io, B.fill)
    print(io, " with indices ", Base.inds2string(inds))
end

# --------------
# ArrayInterface
# --------------

function ArrayInterface.fast_scalar_indexing(A::AlmostBandedMatrix)
    return ArrayInterface.fast_scalar_indexing(typeof(A.bands)) &&
           ArrayInterface.fast_scalar_indexing(typeof(A.fill))
end

function ArrayInterface.qr_instance(A::AlmostBandedMatrix{T}, pivot = NoPivot()) where {T}
    return qr(AlmostBandedMatrix{T}(similar(A.bands, 0, 0), similar(A.fill, 0, 0)))
end

# -------------
# Factorization
# -------------

_factorize(::AbstractAlmostBandedLayout, _, A) = qr(A)
_qr(::AbstractAlmostBandedLayout, _, A) = almostbanded_qr(A)
_qr!(::AlmostBandedLayout, _, A) = almostbanded_qr!(A, Val(false))

almostbanded_qr(A) = _almostbanded_qr(axes(A), A)
function _almostbanded_qr(_, A)
    l, u = almostbandwidths(A)
    B, L = bandpart(A), fillpart(A)
    # Expand the bandsize for the QR factorization
    ## Bypass the safety checks in `AlmostBandedMatrix`
    return almostbanded_qr!(
        AlmostBandedMatrix{eltype(A)}(BandedMatrix(copy(B), (l, l + u)), copy(L)), Val(true))
end

# Band size not yet expanded!
function almostbanded_qr!(R::AbstractMatrix{T}, ::Val{false}) where {T}
    l, u = almostbandwidths(R)
    B, L = bandpart(R), fillpart(R)
    R′ = AlmostBandedMatrix{eltype(R)}(BandedMatrix(B, (l, l + u)), L)
    return almostbanded_qr!(R′, Val(true))
end

function almostbanded_qr!(R::AbstractMatrix{T}, ::Val{true}) where {T}
    return almostbanded_qr!(R, zeros(T, minimum(size(R))))
end

function almostbanded_qr!(R::AbstractMatrix{T}, τ) where {T}
    R′, τ′ = _almostbanded_qr!(R, τ)
    return QR(R′, τ′)
end

function _almostbanded_qr!(A::AbstractMatrix{T}, τ::AbstractVector{T}) where {T}
    m, n = size(A)
    return _almostbanded_qr!(A, τ, min(m - 1 + !(T <: Real), n))
end

@views function _almostbanded_qr!(A::AbstractMatrix, τ::AbstractVector, ncols::Int)
    T = eltype(A)
    B, L = bandpart(A), fillpart(A)
    l, u = bandwidths(B)
    m, n = size(A)
    mf = size(L, 1)
    finish_part_setindex!(B, L)

    U = similar(L, m, mf)
    fill!(U, zero(eltype(L)))
    fill!(U[diagind(U)], one(eltype(L)))

    k = 1
    while k ≤ ncols
        kr = k:min(k + l + u, m)
        jr1 = k:min(k + u, n)
        jr2 = (k + u + 1):min(last(kr) + u, n)
        jr3 = k:min(k + u, n, ncols)
        S = B[kr, jr1]
        τv = τ[jr3]
        R, _ = _banded_qr!(S, τv, length(jr3))
        Q = QRPackedQ(R, τv)

        B_right = B[kr, jr2]
        L_right = L[:, jr2]
        U′ = U[kr, :]
        for j in 1:length(jr2)
            muladd!(
                -one(T), U′[(j + 1):end, :], L_right[:, j], one(T), B_right[(j + 1):end, j])
        end
        banded_qr_lmul!(Q', B_right)
        banded_qr_lmul!(Q', U′)
        for j in 1:length(jr2)
            muladd!(
                one(T), U′[(j + 1):end, :], L_right[:, j], one(T), B_right[(j + 1):end, j])
        end
        k = last(jr1) + 1
    end

    return AlmostBandedMatrix{eltype(A)}(B, Mul(U, L)), τ
end

function getQ(F::QR{<:Any, <:AlmostBandedMatrix})
    return LinearAlgebra.QRPackedQ(bandpart(F.factors), F.τ)
end
function getR(F::QR{<:Any, <:AlmostBandedMatrix})
    n = min(size(F.factors, 1), size(F.factors, 2))
    return UpperTriangular(view(F.factors, 1:n, 1:n))
end

function _almostbanded_longrect_ldiv!(A::QR, B)
    m, n = size(A)
    R = A.factors
    lmul!(adjoint(A.Q), B)
    B̃ = view(B, 1:n, :)
    B̃ .= Ldiv(UpperTriangular(view(R, 1:n, 1:n)), B̃)
    return B
end

function _almostbanded_square_ldiv!(A::QR, B)
    R = A.factors
    lmul!(adjoint(A.Q), B)
    B .= Ldiv(UpperTriangular(R), B)
    return B
end

_almostbanded_widerect_ldiv!(::QR{T}, B) where {T} = error("Not implemented")

const UpperLayoutMatrix{T} = UpperTriangular{T, <:LayoutMatrix{T}}

for Typ in
    (:StridedVector, :StridedMatrix, :AbstractVecOrMat, :UpperLayoutMatrix, :LayoutMatrix)
    @eval function ldiv!(A::QR{T, <:AlmostBandedMatrix}, B::$Typ{T}) where {T}
        m, n = size(A)
        if m == n
            _almostbanded_square_ldiv!(A, B)
        elseif n > m
            _almostbanded_widerect_ldiv!(A, B)
        else
            _almostbanded_longrect_ldiv!(A, B)
        end
    end
end

# needed for adaptive QR
function Base.materialize!(M::Lmul{<:QRPackedQLayout{<:AlmostBandedLayout}})
    return lmul!(QRPackedQ(bandpart(M.A.factors), M.A.τ), M.B)
end
function Base.materialize!(M::Lmul{<:AdjQRPackedQLayout{<:AlmostBandedLayout}})
    Q = M.A'
    return lmul!(QRPackedQ(bandpart(Q.factors), Q.τ)', M.B)
end

triangularlayout(::Type{Tri}, ::ML) where {Tri, ML <: AlmostBandedLayout} = Tri{ML}()

@inline function __arguments(x::LazyArray, ::AlmostBandedMatrix, ::Val)
    return LazyArrays.arguments(x)
end
@inline __arguments(x::Mul, ::AlmostBandedMatrix, ::Val) = (x.A, x.B)
@inline __arguments(x::AbstractArray, ::AlmostBandedMatrix, ::Val{false}) = (nothing, x)
@inline function __arguments(x::AbstractArray, A::AlmostBandedMatrix, ::Val{true})
    U = similar(A, size(A, 1), size(x, 1))
    fill!(U, zero(eltype(A)))
    fill!(U[diagind(U)], one(eltype(A)))
    return U, x
end

@inline __lowrankfillpart(R) = __arguments(fillpart(R), R, Val(false))
@inline @views function __lowrankfillpart(R::SubArray)
    pR = parent(R)
    pU, pV = __lowrankfillpart(pR)
    idx1, idx2 = parentindices(R)
    return pU[idx1, :], pV[:, idx2]
end

@inline __original_almostbandedrank(A) = size(first(__lowrankfillpart(A)), 2)

@views function _almostbanded_upper_ldiv!(
        ::Type{Tri}, R::AbstractMatrix, b::AbstractVector{T}, buffer) where {T, Tri}
    B = bandpart(R)
    U, V = __lowrankfillpart(R)
    fill!(buffer, zero(T))

    l, u = bandwidths(B)
    k = n = size(R, 2)

    while k > 0
        kr = max(1, k - u):k
        jr1 = (k + 1):(k + u + 1)
        jr2 = (k + u + 2):(k + 2u + 2)
        bv = b[kr]
        if jr2[1] < n
            muladd!(one(T), V[:, jr2], b[jr2], one(T), buffer)
            muladd!(-one(T), U[kr, :], buffer, one(T), bv)
        end
        if jr1[1] < n
            muladd!(-one(T), R[kr, jr1], b[jr1], one(T), bv)
        end
        materialize!(Ldiv(Tri(B[kr, kr]), bv))
        k = kr[1] - 1
    end

    return b
end

function Base.materialize!(M::MatLdivVec{TriangularLayout{'U', 'N', AlmostBandedLayout}})
    R, x = M.A, M.B
    A = triangulardata(R)
    r = __original_almostbandedrank(A)
    _almostbanded_upper_ldiv!(UpperTriangular, A, x, Vector{eltype(M)}(undef, r))
    return x
end

function Base.materialize!(M::MatLdivVec{TriangularLayout{'U', 'U', AlmostBandedLayout}})
    R, x = M.A, M.B
    A = triangulardata(R)
    r = __original_almostbandedrank(A)
    _almostbanded_upper_ldiv!(UnitUpperTriangular, A, x, Vector{eltype(M)}(undef, r))
    return x
end

function Base.materialize!(M::MatLdivVec{TriangularLayout{'L', 'N', AlmostBandedLayout}})
    R, x = M.A, M.B
    A = triangulardata(R)
    materialize!(Ldiv(LowerTriangular(bandpart(A)), x))
    return x
end

function Base.materialize!(M::MatLdivVec{TriangularLayout{'L', 'U', AlmostBandedLayout}})
    R, x = M.A, M.B
    A = triangulardata(R)
    materialize!(Ldiv(UnitLowerTriangular(bandpart(A)), x))
    return x
end

# ---------------
# Matrix Multiply
# ---------------

@views function muladd!(
        α, A::AlmostBandedMatrix, B::AbstractVecOrMat, β, C::AbstractVecOrMat)
    L = fillpart(A)
    muladd!(α, L, B, β, selectdim(C, 1, 1:size(L, 1)))
    muladd!(α, exclusive_bandpart(A), B, β, selectdim(C, 1, (size(L, 1) + 1):size(C, 1)))
    return C
end

# ---------------------
# Precompile Some Stuff
# ---------------------

@setup_workload begin
    m = 2
    n = 10
    for T in (Float32, Float64)
        B = brand(T, n, n, m + 1, m)
        F = rand(T, m, n)

        @compile_workload begin
            A = AlmostBandedMatrix(B, F)

            # QR
            fact = qr(A)

            # Linear Solve
            fact \ rand(T, n)
        end
    end
end

export AlmostBandedMatrix, bandpart, fillpart, exclusive_bandpart, finish_part_setindex!,
       almostbandwidths, almostbandedrank

end
