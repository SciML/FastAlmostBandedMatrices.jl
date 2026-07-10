using FastAlmostBandedMatrices
using Test

function _owned_public_api_names()
    source = read(
        joinpath(pkgdir(FastAlmostBandedMatrices), "src", "FastAlmostBandedMatrices.jl"),
        String,
    )
    export_matches = collect(eachmatch(r"(?ms)^export\s+(.+?)(?=\n\n|$)", source))
    exports = Symbol[]
    for export_match in export_matches
        block = replace(export_match.captures[1], "\n" => " ")
        names = filter(!isempty, strip.(split(block, ",")))
        append!(exports, Symbol.(names))
    end
    return exports
end

@testset "public API documentation" begin
    public_names = _owned_public_api_names()
    @test !isempty(public_names)

    missing_docs = Symbol[]
    for name in public_names
        binding = Docs.Binding(FastAlmostBandedMatrices, name)
        Docs.hasdoc(binding) || push!(missing_docs, name)
    end
    @test isempty(missing_docs)

    api_page = read(joinpath(pkgdir(FastAlmostBandedMatrices), "docs", "src", "api.md"), String)
    missing_api_entries = filter(name -> !occursin(String(name), api_page), public_names)
    @test isempty(missing_api_entries)
end
