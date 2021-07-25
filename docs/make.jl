using TopOptLevelSet
using Documenter

DocMeta.setdocmeta!(TopOptLevelSet, :DocTestSetup, :(using TopOptLevelSet); recursive=true)

makedocs(;
    modules=[TopOptLevelSet],
    authors="Mohamed Tarek <mohamed82008@gmail.com> and contributors",
    repo="https://github.com/JuliaTopOpt/TopOptLevelSet.jl/blob/{commit}{path}#{line}",
    sitename="TopOptLevelSet.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaTopOpt.github.io/TopOptLevelSet.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaTopOpt/TopOptLevelSet.jl",
)
