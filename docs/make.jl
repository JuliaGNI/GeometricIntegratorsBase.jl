using GeometricIntegratorsBase
using Documenter

DocMeta.setdocmeta!(GeometricIntegratorsBase, :DocTestSetup, :(using GeometricIntegratorsBase); recursive=true)

makedocs(;
    modules=[GeometricIntegratorsBase],
    authors="Michael Kraus <michael.kraus@ipp.mpg.de> and contributors",
    sitename="GeometricIntegratorsBase.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaGNI.github.io/GeometricIntegratorsBase.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaGNI/GeometricIntegratorsBase.jl",
    devbranch="main",
)
