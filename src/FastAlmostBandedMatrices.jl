module FastAlmostBandedMatrices

import PrecompileTools: @recompile_invalidations, @setup_workload, @compile_workload

@recompile_invalidations begin
    using ArrayInterface, ArrayLayouts, BandedMatrices, ConcreteStructs, LinearAlgebra,
        MatrixFactorizations, Reexport
end

import ArrayLayouts: MemoryLayout, sublayout, sub_materialize, MatLdivVec, materialize!,
    triangularlayout, triangulardata, zero!, _copyto!, colsupport, rowsupport, _qr, _qr!,
    _factorize
import BandedMatrices: _banded_qr!, bandeddata
import LinearAlgebra: ldiv!
import MatrixFactorizations: QR, QRPackedQ, getQ, getR, QRPackedQLayout, AdjQRPackedQLayout

@reexport using BandedMatrices

# ------------------
# AlmostBandedMatrix
# ------------------

abstract type AbstractAlmostBandedLayout <: MemoryLayout end
struct AlmostBandedLayout <: AbstractAlmostBandedLayout end

@concrete struct AlmostBandedMatrix{T} <: LayoutMatrix{T}
    bands
    fill
end

function AlmostBandedMatrix(fill::AbstractMatrix, bands::BandedMatrix)
    @assert size(fill, 2) == size(bands, 2)
    T = promote_type(eltype(fill), eltype(bands))
    @assert bandwidths(bands)[1] ≥ size(fill, 1) - 1
    @assert ArrayInterface.can_setindex(bands)
    @assert ArrayInterface.fast_scalar_indexing(fill)&&ArrayInterface.fast_scalar_indexing(bands) "Non Fast Scalar Index-able Arrays are not supported currently."
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

@inline bandpart(A::AlmostBandedMatrix) = A.bands
@inline fillpart(A::AlmostBandedMatrix) = A.fill

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

@inline function colsupport(::AbstractAlmostBandedLayout, A, j)
    l, _ = almostbandwidths(A)
    return Base.OneTo(min(maximum(j) + l, size(A, 1)))
end

@inline function rowsupport(::AbstractAlmostBandedLayout, A, k)
    l, _ = almostbandwidths(A)
    return max(1, minimum(k) - l):size(A, 2)
end

@inline function Base.getindex(B::AlmostBandedMatrix, k::Integer, j::Integer)
    if j > k + bandwidth(B.bands, 2) && k ≤ size(B.fill, 1)
        return @inbounds B.fill[k, j]
    else
        return @inbounds B.bands[k, j]
    end
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

# --------------
# ArrayInterface
# --------------

function ArrayInterface.fast_scalar_indexing(A::AlmostBandedMatrix)
    return ArrayInterface.fast_scalar_indexing(typeof(A.bands)) &&
           ArrayInterface.fast_scalar_indexing(typeof(A.fill))
end

function ArrayInterface.qr_instance(A::AlmostBandedMatrix{T}, pivot = NoPivot()) where {T}
    return qr(AlmostBandedMatrix{T}(brand(T, 2, 2, 1, 1), similar(A.fill, 1, 2)))
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
    return almostbanded_qr!(AlmostBandedMatrix{eltype(A)}(BandedMatrix(copy(B), (l, l + u)),
            copy(L)), Val(true))
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

# TODO: We can optimize the memory usage further, but we can do it once it becomes a
#       bottleneck
@views function _almostbanded_qr!(A::AbstractMatrix, τ::AbstractVector, ncols::Int)
    B, L_ = bandpart(A), fillpart(A)
    l, u = bandwidths(B)
    m, n = size(A)
    mf = size(L_, 1)
    finish_part_setindex!(B, L_)

    # Materialize the fill part into a complete matrix
    L = similar(L_, m, n)
    copyto!(L[1:size(L_, 1), :], L_)
    fill!(L[(size(L_, 1) + 1):end, :], zero(eltype(L_)))

    ## L = U * L_
    U = similar(L_, m, mf)
    fill!(U, zero(eltype(L_)))
    fill!(U[diagind(U)], one(eltype(L_)))

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
        L_right = L[kr, jr2]
        for j in 1:length(jr2)
            B_right[(j + 1):end, j] .-= L_right[(j + 1):end, j]
        end
        lmul!(Q', B_right)
        lmul!(Q', U[kr, :])
        mul!(L_right, U[kr, :], L_[:, jr2])
        for j in 1:length(jr2)
            B_right[(j + 1):end, j] .+= L_right[(j + 1):end, j]
        end
        k = last(jr1) + 1
    end

    return AlmostBandedMatrix{eltype(A)}(B, L[1:(m - u), :]), τ
end

function getQ(F::QR{<:Any, <:AlmostBandedMatrix})
    return LinearAlgebra.QRPackedQ(bandpart(F.factors), F.τ)
end
function getR(F::QR{<:Any, <:AlmostBandedMatrix})
    n = min(size(F.factors, 1), size(F.factors, 2))
    return UpperTriangular(view(F.factors, 1:n, 1:n))
end

function ldiv!(A::QR{<:Any, <:AlmostBandedMatrix}, B::AbstractVecOrMat)
    R = A.factors
    lmul!(adjoint(A.Q), B)
    B .= Ldiv(UpperTriangular(R), B)
    return B
end

# needed for adaptive QR
function Base.materialize!(M::Lmul{<:QRPackedQLayout{AlmostBandedLayout}})
    return lmul!(QRPackedQ(bandpart(M.A.factors), M.A.τ), M.B)
end
function Base.materialize!(M::Lmul{<:AdjQRPackedQLayout{AlmostBandedLayout}})
    Q = M.A'
    return lmul!(QRPackedQ(bandpart(Q.factors), Q.τ)', M.B)
end

# TODO: Upper Triangular `ldiv!`
triangularlayout(::Type{Tri}, ::ML) where {Tri, ML <: AlmostBandedLayout} = Tri{ML}()

# -------------
# LinearAlgebra
# -------------

export AlmostBandedMatrix, bandpart, fillpart, exclusive_bandpart, finish_part_setindex!

end
