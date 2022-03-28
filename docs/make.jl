using RandomStrainDistributions
using Documenter

DocMeta.setdocmeta!(RandomStrainDistributions, :DocTestSetup, :(using RandomStrainDistributions); recursive=true)

makedocs(;
    modules=[RandomStrainDistributions],
    authors="W. Joe Meese <meese022@umn.edu> and contributors",
    repo="https://github.com/meesewj/RandomStrainDistributions.jl/blob/{commit}{path}#{line}",
    sitename="RandomStrainDistributions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
