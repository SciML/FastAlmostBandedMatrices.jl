using SciMLTesting, FastAlmostBandedMatrices

const REEXPORTED_API = (
    :Band, :BandError, :BandRange, :BandedMatrices, :BandedMatrix,
    :Eye, :Fill, :Ones, :Zeros,
    :band, :bandrange, :bandwidth, :bandwidths, :brand, :brandn,
    :colrange, :rowrange, :symrcm,
)

run_qa(
    FastAlmostBandedMatrices;
    explicit_imports = true,
    api_docs_kwargs = (;
        rendered = true,
        ignore = REEXPORTED_API,
        rendered_ignore = REEXPORTED_API,
    ),
    # 19 method ambiguities, all in FastAlmostBandedMatrices' own ldiv!/__arguments
    # against ArrayLayouts/LinearAlgebra Triangular/Factorization methods.
    # https://github.com/SciML/FastAlmostBandedMatrices.jl/issues/71
    aqua_broken = (:ambiguities,),
    ei_kwargs = (;
        # Non-public names this package legitimately extends/uses from upstream:
        # ArrayLayouts MatLdivVec/sublayout/triangulardata/triangularlayout/_qr/_qr!/
        #   _factorize/QRPackedQLayout/AdjQRPackedQLayout, MatrixFactorizations
        #   QR/QRPackedQ/getQ/getR, BandedMatrices _banded_qr!/banded_qr_lmul!.
        all_explicit_imports_are_public = (;
            ignore = (
                :MatLdivVec, :sublayout, :triangulardata, :triangularlayout,
                :_qr, :_qr!, :_factorize, :QRPackedQLayout, :AdjQRPackedQLayout,
                :QR, :QRPackedQ, :getQ, :getR,
                :_banded_qr!, :banded_qr_lmul!,
            ),
        ),
        # Qualified accesses of non-public names: Base OneTo/array_summary/dims2string/
        #   inds2string/materialize!, LinearAlgebra QRPackedQ, LazyArrays arguments,
        #   ArrayInterface fast_scalar_indexing/qr_instance.
        all_qualified_accesses_are_public = (;
            ignore = (
                :OneTo, :array_summary, :dims2string, :inds2string, :materialize!,
                :QRPackedQ, :arguments, :fast_scalar_indexing, :qr_instance,
            ),
        ),
    ),
    # 31 names implicitly imported via the package's `using ArrayInterface, ArrayLayouts,
    # BandedMatrices, ConcreteStructs, LazyArrays, LinearAlgebra, ...` plus
    # `@reexport using BandedMatrices`; explicit-import conversion tracked separately.
    # https://github.com/SciML/FastAlmostBandedMatrices.jl/issues/71
    ei_broken = (:no_implicit_imports,),
)
