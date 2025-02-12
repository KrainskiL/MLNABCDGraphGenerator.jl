
struct ABCDParams
    w::Vector{Int}
    s::Vector{Int}
    ξ::Float64
    μ::Nothing
    isCL::Bool
    islocal::Bool
    hasoutliers::Bool

    function ABCDParams(w::Vector{Int}, s::Vector{Int}, ξ::Float64)
        length(w) == sum(s) || throw(ArgumentError("inconsistent data"))
        0 ≤ ξ ≤ 1 || throw(ArgumentError("inconsistent data ξ"))

        new(w, s, ξ, nothing, false, false, false)
    end
end

mutable struct MLNConfig
    seed::Union{Int,Nothing}
    n::Int
    edges_cor::String
    layer_params::String
    d_max_iter::Int
    c_max_iter::Int
    t::Int
    eps::Float64
    d::Int
    edges_filename::String
    communities_filename::String
    l::Int
    qs::Vector{Float64}
    ns::Vector{Int}
    taus::Vector{Float64}
    rs::Vector{Float64}
    gammas::Vector{Float64}
    d_mins::Vector{Int}
    d_maxs::Vector{Int}
    betas::Vector{Float64}
    c_mins::Vector{Int}
    c_maxs::Vector{Int}
    xis::Vector{Float64}
    skip_edges_correlation::Bool
    edges_cor_matrix::Matrix{Float64}
end

"""
    MLNConfig(
        seed::Union{Int,Nothing},
        n::Int,
        edges_cor::String,
        layer_params::String,
        d_max_iter::Int,
        c_max_iter::Int,
        t::Int,
        eps::Float64,
        d::Int,
        edges_filename::String,
        communities_filename::String)
        ::MLNConfig

Creates MLNABCD configuration struct.

**Arguments**
* `filename::String` configuration file name
"""
function MLNConfig(
    seed::Union{Int,Nothing},
    n::Int,
    edges_cor::String,
    layer_params::String,
    d_max_iter::Int,
    c_max_iter::Int,
    t::Int,
    eps::Float64,
    d::Int,
    edges_filename::String,
    communities_filename::String)::MLNConfig

    isnothing(seed) || Random.seed!(seed)
    lparams = readdlm(layer_params, ',', Float64, skipstart=1)
    l = size(lparams)[1]
    skip_edges_correlation = false
    if isempty(edges_cor)
        skip_edges_correlation = true
        edges_cor_matrix = zeros(Float64, l, l)
    elseif occursin(".csv", edges_cor)
        edges_cor_matrix = convert.(Float64, readdlm(edges_cor, ',', skipstart=1)[:, 2:end])
        x, y = size(edges_cor_matrix)
        @assert (x == l) && (y == l) "Edges correlation matrix size don't match number of layers from $(layer_params)"
    else
        edges_cor_matrix = fill(parse(Float64, edges_cor), l, l)
    end
    qs = lparams[:, 1]
    MLNConfig(seed, n, edges_cor, layer_params, d_max_iter, c_max_iter,
        t, eps, d, edges_filename, communities_filename,
        l, qs, round.(Int, n .* qs), lparams[:, 2], lparams[:, 3],
        lparams[:, 4], Int.(lparams[:, 5]), Int.(lparams[:, 6]), lparams[:, 7],
        Int.(lparams[:, 8]), Int.(lparams[:, 9]), lparams[:, 10],
        skip_edges_correlation, edges_cor_matrix)
end

"""
    parse_config(filename::String)::MLNConfig

Parse configuration file into MLNABCD configuration struct.

**Arguments**
* `filename::String` configuration file name
"""
function parse_config(filename::String)::MLNConfig
    conf = Pkg.TOML.parsefile(filename)
    return MLNConfig(
        length(conf["seed"]) > 0 ? parse(Int, conf["seed"]) : nothing,
        parse(Int, conf["n"]),
        conf["edges_cor"],
        conf["layer_params"],
        parse(Int, conf["d_max_iter"]),
        parse(Int, conf["c_max_iter"]),
        parse(Int, conf["t"]),
        parse(Float64, conf["e"]),
        parse(Int, conf["d"]),
        haskey(conf, "edges_filename") ? conf["edges_filename"] : "edges.dat",
        haskey(conf, "communities_filename") ? conf["communities_filename"] : "communities.dat")
end

"""
    generate_active_nodes(cfg::MLNConfig)::Vector{Vector{Int}}

Samples active nodes for each layer.

**Arguments**
* `cfg::MLNConfig` MLNABCD configuration
"""
function generate_active_nodes(cfg::MLNConfig)::Vector{Vector{Int}}
    return sample.(Ref(1:cfg.n), cfg.ns, replace=false, ordered=true)
end

"""
    generate_abcd(cfg::MLNConfig, degs::Vector{Vector{Int}},
        coms_sizes::Vector{Vector{Int}}, coms::Vector{Vector{Int}})
        ::Vector{Set{Tuple{Int,Int}}}

Generates edges of multilayer network using [ABCD](https://github.com/bkamins/ABCDGraphGenerator.jl) model.

**Arguments**
* `cfg::MLNConfig` MLNABCD configuration
* `degs::Vector{Vector{Int}}` degree sequences for active nodes
* `coms_sizes::Vector{Vector{Int}}` community sizes
* `coms::Vector{Vector{Int}}` assignment of active nodes to communities
"""
function generate_abcd(cfg::MLNConfig,
    degs::Vector{Vector{Int}},
    coms_sizes::Vector{Vector{Int}},
    coms::Vector{Vector{Int}})::Vector{Set{Tuple{Int,Int}}}
    p = MLNABCDGraphGenerator.ABCDParams.(degs, coms_sizes, cfg.xis)
    return ABCDGraphGenerator.config_model.(coms, p)
end

"""
    write_edges(cfg::MLNConfig, edges::Vector{Set{Tuple{Int,Int}}})

Saves edges of multilayer network in `cfg.edges_filename` file.

**Arguments**
* `cfg::MLNConfig` MLNABCD configuration
* `edges::Vector{Set{Tuple{Int,Int}}}` edges between agents in each layer
"""
function write_edges(cfg::MLNConfig, edges::Vector{Set{Tuple{Int,Int}}})
    open(cfg.edges_filename, "w") do io
        for i in eachindex(edges)
            for (a, b) in sort!(collect(edges[i]))
                println(io, a, "\t", b, "\t", i)
            end
        end
    end
end

"""
    write_communities(cfg::MLNConfig, coms::Vector{Vector{Int}})

Saves communities of multilayer network in `cfg.communities_filename` file.

**Arguments**
* `cfg::MLNConfig` MLNABCD configuration
* `coms::Vector{Vector{Int}}` assignment of agents to communities for each layer
"""
function write_communities(cfg::MLNConfig, coms::Vector{Vector{Int}})
    open(cfg.communities_filename, "w") do io
        for i in eachindex(coms)
            for c in coms[i]
                println(io, c, "\t", i)
            end
        end
    end
end