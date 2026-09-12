using SciMLTesting, FastAlmostBandedMatrices

# The BandedMatrices.jl names FastAlmostBandedMatrices intentionally reexports, because they
# are needed for its own documented use: `AlmostBandedMatrix(bands::BandedMatrix, fill)` is
# the documented constructor and the documented examples build `bands` with `brand`.
# Kept in sync with the reexport `export` block in src/FastAlmostBandedMatrices.jl.
const REEXPORTED_API = (
    :Band,
    :BandError,
    :BandRange,
    :BandedMatrix,
    :band,
    :bandrange,
    :bandwidth,
    :bandwidths,
    :brand,
    :brandn,
    :colrange,
    :rowrange,
)

run_qa(
    FastAlmostBandedMatrices;
    reexports_allow = REEXPORTED_API,
    ei_kwargs = (;
        # Non-public names this package legitimately extends/uses from upstream:
        # ArrayLayouts MatLdivVec/sublayout/triangulardata/triangularlayout/_qr/_qr!/
        #   _factorize/QRPackedQLayout/AdjQRPackedQLayout, MatrixFactorizations
        #   QR/QRPackedQ/getQ/getR, BandedMatrices _banded_qr!/banded_qr_lmul!.
        all_explicit_imports_are_public = (;
            ignore = (
                :MatLdivVec,
                :sublayout,
                :triangulardata,
                :triangularlayout,
                :_qr,
                :_qr!,
                :_factorize,
                :QRPackedQLayout,
                :AdjQRPackedQLayout,
                :QR,
                :QRPackedQ,
                :getQ,
                :getR,
                :_banded_qr!,
                :banded_qr_lmul!,
            ),
        ),
        # Qualified accesses of non-public names: Base OneTo/array_summary/dims2string/
        #   inds2string/materialize!, LinearAlgebra QRPackedQ, LazyArrays arguments,
        #   ArrayInterface fast_scalar_indexing/qr_instance.
        all_qualified_accesses_are_public = (;
            ignore = (
                :OneTo,
                :array_summary,
                :dims2string,
                :inds2string,
                :materialize!,
                :QRPackedQ,
                :arguments,
                :fast_scalar_indexing,
                :qr_instance,
            ),
        ),
    ),
)
