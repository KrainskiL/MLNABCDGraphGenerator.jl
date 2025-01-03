using ABCDGraphGenerator
using DelimitedFiles
using MLNABCDGraphGenerator
using Pkg
using Random
using StatsBase

@info "Usage: julia mlnabcd_sampler.jl config_filename"
@info "For the syntax of config_filename see example_config.toml file"

filename = ARGS[1]
conf = Pkg.TOML.parsefile(filename)
isempty(conf["seed"]) || Random.seed!(parse(Int, conf["seed"]))

n = parse(Int, conf["n"])
t = parse(Int, conf["t"])
ϵ = parse(Float64, conf["e"])
islocal = false
isCL = false
μ = nothing

lparams_file = conf["layer_params"]
lparams = readdlm(lparams_file, ',', skipstart=1)
l = size(lparams)[1]
qs = lparams[:, 1]
taus = lparams[:, 2]
rs = lparams[:, 3]

skip_edges_correlation = false
if isempty(conf["edges_cor"]) || conf["edges_cor"] == ""
    skip_edges_correlation = true
elseif occursin(".csv", conf["edges_cor"])
    edges_cor_matrix = convert.(Float64,readdlm(conf["edges_cor"],',', skipstart=1)[:,2:end])
    x, y = size(edges_cor_matrix)
    @assert (x == l) && (y == l) "Edges correlation matrix size don't match number of layers from $(lparams_file)"
else
    edges_cor = parse(Float64, conf["edges_cor"])
    edges_cor_matrix = fill(edges_cor,l,l)
end

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
degs = MLNABCDGraphGenerator.degrees_correlation(n, degs, taus, active_nodes, false)

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
if !skip_edges_correlation
    MLNABCDGraphGenerator.adjust_edges_correlation!(edges, coms, active_nodes, edges_cor_matrix, t, ϵ, false, false)
end
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