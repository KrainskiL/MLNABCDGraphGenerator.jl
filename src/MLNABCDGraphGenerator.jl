module MLNABCDGraphGenerator

using Random
using StatsBase
using ArgParse
using ABCDGraphGenerator
using DelimitedFiles

include("degree_cor.jl")
include("community_cor.jl")
include("edges_cor.jl")

end # module MLNABCDGraphGenerator