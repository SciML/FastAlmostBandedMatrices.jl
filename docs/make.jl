using Documenter, FastAlmostBandedMatrices

makedocs(;
    sitename = "FastAlmostBandedMatrices.jl",
    authors = "Avik Pal et al.",
    modules = [FastAlmostBandedMatrices],
    clean = true,
    doctest = false,
    linkcheck = false,
    checkdocs = :exports,
    format = Documenter.HTML(;
        canonical = "https://docs.sciml.ai/FastAlmostBandedMatrices/stable/"
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ]
)

deploydocs(;
    repo = "github.com/SciML/FastAlmostBandedMatrices.jl.git",
    push_preview = true
)
