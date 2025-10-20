using Documenter
using DocumenterInterLinks
using GeometricIntegratorsBase
using GeometricEquations
using GeometricSolutions

DocMeta.setdocmeta!(GeometricIntegratorsBase, :DocTestSetup, :(using GeometricIntegratorsBase); recursive=true)

links = InterLinks(
    "GeometricEquations" => "https://JuliaGNI.github.io/GeometricEquations.jl/stable/",
)

makedocs(;
    sitename="GeometricIntegratorsBase.jl",
    plugins=[links,],
    warnonly=Documenter.except(:autodocs_block, :cross_references, :docs_block, :doctest, :eval_block, :example_block, :footnote, :linkcheck_remotes, :linkcheck, :meta_block, :parse_error, :setup_block),
    authors="Michael Kraus <michael.kraus@ipp.mpg.de> and contributors",
    format=Documenter.HTML(;
        canonical="https://JuliaGNI.github.io/GeometricIntegratorsBase.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Extrapolation Methods" => "extrapolation.md",
        "Euler Methods" => "euler.md",
        "Dependencies" => [
            "Equations" => "deps/equations.md",
            "Problems" => "deps/problems.md",
            "Solutions" => "deps/solutions.md",
        ],
    ],
    modules=[
        GeometricIntegratorsBase,
        GeometricEquations,
        GeometricSolutions,
    ],
)

deploydocs(;
    repo="github.com/JuliaGNI/GeometricIntegratorsBase.jl",
    devbranch="main",
)
