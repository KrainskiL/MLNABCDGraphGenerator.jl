using ABCDGraphGenerator
using DelimitedFiles
using MLNABCDGraphGenerator
using Pkg
using Random
using StatsBase

@info "Usage: julia mlnabcd_sampler.jl config_filename"
@info "For the syntax of config_filename see example_config.toml file"

filename = ARGS[1]
filename = "example_config.toml"
conf = Pkg.TOML.parsefile(filename)
isempty(conf["seed"]) || Random.seed!(parse(Int, conf["seed"]))

n = parse(Int, conf["n"])
t = parse(Int, conf["t"])
ϵ = parse(Float64, conf["e"])
edges_r = parse(Float64, conf["R"])
islocal = false
isCL = false
μ = nothing

lparams = readdlm(conf["layer_params"], ',', skipstart=1)
l = size(lparams)[1]
qs = lparams[:, 1]
rhos = lparams[:, 2]
rs = lparams[:, 3]

#Active nodes
ns = round.(Int, n .* qs)
active_nodes = [sample(1:n, ni, replace=false, ordered=true) for ni in ns]

#Degree Sequences
d_max_iter = parse(Int, conf["d_max_iter"])
degs = Vector{Vector{Int}}()
for i in 1:size(lparams)[1]
    τ, d_min, d_max = lparams[i, 4:6]
    ni = ns[i]
    degs_layer = ABCDGraphGenerator.sample_degrees(τ, Int(d_min), Int(d_max), ni, d_max_iter)
    push!(degs, degs_layer)
end
degs = MLNABCDGraphGenerator.degrees_correlation(n, degs, rhos, active_nodes, false)
# for i in eachindex(degs_correlated)
#     open(io -> foreach(d -> println(io, d), degs_correlated[i]), "deg_$(i).dat", "w")
# end

#Sizes of communities
c_max_iter = parse(Int, conf["c_max_iter"])
coms_sizes = []
for i in 1:size(lparams)[1]
    τ, c_min, c_max = lparams[i, 7:9]
    ni = ns[i]
    coms_layer = ABCDGraphGenerator.sample_communities(τ, Int(c_min), Int(c_max), ni, c_max_iter)
    push!(coms_sizes, coms_layer)
end
coms = MLNABCDGraphGenerator.communities_correlation(n, coms_sizes, rs, active_nodes)
# for i in eachindex(coms)
#     open(io -> foreach(d -> println(io, d), coms[i]), "comm_$(i).dat", "w")
# end

edges = Vector{Set{Tuple{Int64,Int64}}}()
clusters = []
for i in 1:l
    ξ = lparams[i, 10]
    p = MLNABCDGraphGenerator.MLNABCDParams(degs[i], coms_sizes[i], ξ)
    edges_layer = ABCDGraphGenerator.config_model(coms[i], p)
    push!(edges, edges_layer)
end
edges = MLNABCDGraphGenerator.map_edges_to_agents(edges, active_nodes)
coms = MLNABCDGraphGenerator.map_communities_to_agents(n, coms, active_nodes)
MLNABCDGraphGenerator.adjust_edges_correlation!(edges, coms, active_nodes, edges_r, t, ϵ, false)
edges_common_agents = MLNABCDGraphGenerator.common_agents_edges(edges, active_nodes)
edges_cor = MLNABCDGraphGenerator.calculate_edges_cor(edges_common_agents, true)

open("edges.dat", "w") do io
    for i in eachindex(edges)
        for (a, b) in sort!(collect(edges[i]))
            println(io, a, "\t", b, "\t", i)
        end
    end
end
open("communities.dat", "w") do io
    for i in eachindex(coms)
        for c in coms[i]
            println(io, c, "\t", i)
        end
    end
end