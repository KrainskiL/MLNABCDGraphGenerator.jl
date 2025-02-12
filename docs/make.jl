using Documenter

try
    using MLNABCDGraphGenerator
catch
    if !("../src/" in LOAD_PATH)
       push!(LOAD_PATH,"../src/")
       @info "Added \"../src/\"to the path: $LOAD_PATH "
       using MLNABCDGraphGenerator
    end
end

makedocs(
    sitename = "MLNABCDGraphGenerator.jl",
    format = format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    modules = [MLNABCDGraphGenerator],
    pages = ["index.md", "reference.md"],
    doctest = false
)



deploydocs(
    repo ="github.com/KrainskiL/MLNABCDGraphGenerator.jl.git",
    target="build"
)