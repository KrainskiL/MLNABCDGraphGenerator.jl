using Documenter, MLNABCDGraphGenerator

makedocs(
    sitename="MLNABCDGraphGenerator.jl",
    format=format = Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true"
    ),
    modules=[MLNABCDGraphGenerator],
    pages=["index.md", "reference.md"],
    doctest=false
)

deploydocs(
    repo="github.com/KrainskiL/MLNABCDGraphGenerator.jl.git",
    target="build"
)