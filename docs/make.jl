using CombinedUncertainDiffEq
using Documenter

DocMeta.setdocmeta!(CombinedUncertainDiffEq, :DocTestSetup, :(using CombinedUncertainDiffEq); recursive=true)

makedocs(;
    modules=[CombinedUncertainDiffEq],
    authors="Jonnie Diegelman <jonniediegelman@gmail.com> and contributors",
    repo="https://github.com/jonniedie/CombinedUncertainDiffEq.jl/blob/{commit}{path}#{line}",
    sitename="CombinedUncertainDiffEq.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jonniedie.github.io/CombinedUncertainDiffEq.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jonniedie/CombinedUncertainDiffEq.jl",
)
