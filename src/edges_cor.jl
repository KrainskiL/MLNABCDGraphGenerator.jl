
"""
    map_edges_to_agents(edges::Vector{Set{Tuple{Int64,Int64}}}, active_nodes::Vector{Vector{Int64}})

Translates edges to agent-based representation.

**Arguments**
* `edges::Vector{Set{Tuple{Int64,Int64}}}` edges based on layer-specific node labels
* `active_nodes::Vector{Vector{Int64}}` IDs of active agents
"""
function map_edges_to_agents(edges::Vector{Set{Tuple{Int64,Int64}}}, active_nodes::Vector{Vector{Int64}})
    edges_agents = Vector{Set{Tuple{Int64,Int64}}}()
    for (i, layer) in enumerate(edges)
        push!(edges_agents, Set([(active_nodes[i][edge[1]], active_nodes[i][edge[2]]) for edge in layer]))
    end
    return edges_agents
end

"""
    map_communities_to_agents(
        n_agents::Int,
        communities::Vector{Vector{Int64}},
        active_nodes::Vector{Vector{Int64}})

Translates communities assignment to agent-based representation and imputes 0 for non-active agents.

**Arguments**
* `n_agents::Int` number of agents
* `communities::Vector{Vector{Int64}}` assignment of nodes to communities
* `active_nodes::Vector{Vector{Int64}}` IDs of active agents
"""
function map_communities_to_agents(n_agents::Int, communities::Vector{Vector{Int64}}, active_nodes::Vector{Vector{Int64}})
    l = length(communities)
    coms_agents = [zeros(Int, n_agents) for i in 1:l]
    for i in 1:l
        for j in eachindex(communities[i])
            coms_agents[i][active_nodes[i][j]] = communities[i][j]
        end
    end
    return coms_agents
end

"""
    common_agents_dict(active_nodes::Vector{Vector{Int64}})

Finds common active agents between all pairs of layers.

**Arguments**
* `active_nodes::Vector{Vector{Int64}}` IDs of active agents
"""
function common_agents_dict(active_nodes::Vector{Vector{Int64}})
    common_agents_dict = Dict{Tuple{Int,Int},Set{Int}}()
    for i in eachindex(active_nodes)
        for j in (i+1):length(active_nodes)
            common_agents_dict[(i, j)] = BitSet(intersect(Set(active_nodes[i]), Set(active_nodes[j])))
        end
    end
    return common_agents_dict
end

"""
    common_agents_edges(
        edges_agents::Vector{Set{Tuple{Int64,Int64}}},
        common_agents_dict::Dict{Tuple{Int,Int},Set{Int}})

Finds edges including only common active agents for all pairs of layers.

**Arguments**
* `edges_agents::Vector{Set{Tuple{Int64,Int64}}}` agent-based edges of multilayer network
* `common_agents_dict::Dict{Tuple{Int,Int},Set{Int}}` mapping of common active agents between layers
"""
function common_agents_edges(edges_agents::Vector{Set{Tuple{Int64,Int64}}}, common_agents_dict::Dict{Tuple{Int,Int},Set{Int}})
    edges_on_common_agents_dict = Dict{Tuple{Int,Int},Vector{Set{Tuple{Int,Int}}}}()
    for i in eachindex(edges_agents)
        for j in (i+1):length(edges_agents)
            common_agents = common_agents_dict[(i, j)]
            edges_i = Set([e for e in edges_agents[i] if issubset(e, common_agents)])
            edges_j = Set([e for e in edges_agents[j] if issubset(e, common_agents)])
            edges_on_common_agents_dict[(i, j)] = [copy(edges_i), copy(edges_j)]
        end
    end
    return edges_on_common_agents_dict
end

"""
    map_neighbours(edges::Set{Tuple{Int,Int}}, coms::Vector{Int})

Finds neigbhours of each agent within and outside of given node community.

**Arguments**
* `edges::Set{Tuple{Int,Int}}` agent-based edges of multilayer network
* `coms::Vector{Int}` assignment of nodes to communities
"""
function map_neighbours(edges::Set{Tuple{Int,Int}}, coms::Vector{Int})
    nmap = Dict{Int,Set{Int}}()
    for (src, dst) in edges
        mult = coms[src] == coms[dst] ? 1 : -1 #positive for communities, negative for background
        for (key, val) in [(mult * src, dst), (mult * dst, src)]
            if haskey(nmap, key)
                push!(nmap[key], val)
            else
                nmap[key] = Set(val)
            end
        end
    end
    return nmap
end

"""
    map_communities(edges::Set{Tuple{Int,Int}}, coms::Vector{Int})

Finds edges within communities and background graph.

**Arguments**
* `edges::Set{Tuple{Int,Int}}` agent-based edges of multilayer network
* `coms::Vector{Int}` assignment of nodes to communities
"""
function map_communities(edges::Set{Tuple{Int,Int}}, coms::Vector{Int})
    nmap = Dict{Int,Set{Tuple{Int,Int}}}()
    for (src, dst) in edges
        src_com = coms[src]
        dst_com = coms[dst]
        key = src_com == dst_com ? src_com : -1
        if haskey(nmap, key)
            push!(nmap[key], (src, dst))
        else
            nmap[key] = Set([(src, dst)])
        end
    end
    return nmap
end

"""
    neighbours_common_agents_edges(
        edges_common_agents::Dict{Tuple{Int,Int},Vector{Set{Tuple{Int,Int}}}},
        coms_agents::Vector{Vector{Int64}})

Finds neighbours of each agent in graphs induced by common agents between layers.

**Arguments**
* `edges_common_agents::Dict{Tuple{Int,Int},Vector{Set{Tuple{Int,Int}}}}` graphs induced by common agents between layers
* `coms_agents::Vector{Vector{Int64}}` assignment of nodes to communities in each layer
"""
function neighbours_common_agents_edges(
    edges_common_agents::Dict{Tuple{Int,Int},Vector{Set{Tuple{Int,Int}}}},
    coms_agents::Vector{Vector{Int64}})

    nbhs_edges_on_common_agents_dict = Dict()
    for (k, v) in edges_common_agents
        nbhs_edges_on_common_agents_dict[k] = (map_neighbours(v[1], coms_agents[k[1]]), map_neighbours(v[2], coms_agents[k[2]]))
    end
    return nbhs_edges_on_common_agents_dict
end

"""
    edges_in_communities_common_agents_edges(
        edges_common_agents::Dict{Tuple{Int,Int},Vector{Set{Tuple{Int,Int}}}},
        coms_agents::Vector{Vector{Int64}})

Finds edges within communities and background graph in graphs induced by common agents between layers.

**Arguments**
* `edges_common_agents::Dict{Tuple{Int,Int},Vector{Set{Tuple{Int,Int}}}}` graphs induced by common agents between layers
* `coms_agents::Vector{Vector{Int64}}` assignment of nodes to communities in each layer
"""
function edges_in_communities_common_agents_edges(
    edges_common_agents::Dict{Tuple{Int,Int},Vector{Set{Tuple{Int,Int}}}},
    coms_agents::Vector{Vector{Int64}})

    edges_in_coms = Dict()
    for (k, v) in edges_common_agents
        edges_in_coms[k] = (map_communities(v[1], coms_agents[k[1]]), map_communities(v[2], coms_agents[k[2]]))
    end
    return edges_in_coms
end

"""
    calculate_edges_cor(
        edges_common_agents::Dict{Tuple{Int,Int},Vector{Set{Tuple{Int,Int}}}},
        desired_cor::Matrix{Float64})

Calculates correlation between edges in each pair of layers, difference between desired correlation and distance between desired correlation matrix and current correlation matrix.

**Arguments**
* `edges_common_agents::Dict{Tuple{Int,Int},Vector{Set{Tuple{Int,Int}}}}` graphs induced by common agents between layers
* `desired_cor::Matrix{Float64}` input correlation matrix between layers
"""
function calculate_edges_cor(
    edges_common_agents::Dict{Tuple{Int,Int},Vector{Set{Tuple{Int,Int}}}},
    desired_cor::Matrix{Float64})

    l = size(desired_cor)[1]
    current_rs = deepcopy(desired_cor)
    for (k, v) in edges_common_agents
        e1, e2 = v
        r = length(intersect(e1, e2)) / (min(length(e1), length(e2)))
        current_rs[k...] = isnan(r) ? 0.0 : r
    end
    rs_diff = current_rs - desired_cor
    rs_dist = 0
    for i in 1:l
        for j in (i+1):l
            rs_dist += rs_diff[i, j]^2
        end
    end
    return (current_rs, rs_diff, rs_dist)
end

"""
    pick_cor_to_improve(cor_diff::Matrix{Float64}, method::String="weighted")

Selects entry of correlation matrix to improve in next batch.

**Arguments**
* `cor_diff::Matrix{Float64}` difference between desired and current correlation matrices
* `method::String="weighted"` method of entry selection
"""
function pick_cor_to_improve(cor_diff::Matrix{Float64}, method::String="weighted")
    @assert method in ["weighted", "uniform", "e-greedy"] "Method must be weighted, uniform or e-greedy, got $(method)"
    l = size(cor_diff)[1]
    cor_diff_vec = []
    for i in 1:l
        for j in (i+1):l
            push!(cor_diff_vec, ((i, j), cor_diff[i, j]))
        end
    end
    cor_diff_values = abs.(getindex.(cor_diff_vec, 2))
    if method == "uniform"
        pick = rand(cor_diff_vec)
    elseif method == "weighted"
        pick = sample(cor_diff_vec, Weights(cor_diff_values))
    else
        ϵ = 0.5
        max_diff_idx = argmax(cor_diff_values)
        other_weight = (1 - ϵ) / (length(cor_diff_vec) - 1)
        weights = fill(other_weight, length(cor_diff_vec))
        weights[max_diff_idx] = ϵ
        pick = sample(cor_diff_vec, Weights(weights))
    end
    layers, diff = pick
    if diff > 0
        direction = "decrease"
    elseif diff < 0
        direction = "increase"
    else
        direction = "skip"
    end
    return (layers, direction)
end

"""
    write_cor_diff_to_file(io::IOStream, cor_diff::Matrix{Float64}, cor_dist::Float64)

Writes differences between current and desired correlation matrices entries to file. Used to track changes in empirical correlation between batches.

**Arguments**
* `io::IOStream` IO handler to a file
* `cor_diff::Matrix{Float64}` difference between desired and current correlation matrices
* `cor_dist::Float64` distance between desired and current correlation matrices
"""
function write_cor_diff_to_file(io::IOStream, cor_diff::Matrix{Float64}, cor_dist::Float64)
    l = size(cor_diff)[1]
    line = ""
    for i in 1:l
        for j in (i+1):l
            line *= string(cor_diff[i, j]) * ","
        end
    end
    println(io, line * string(cor_dist))
end

"""
    increase_edges_correlation!(layers::Tuple{Int,Int},
        edges::Vector{Set{Tuple{Int,Int}}},
        edges_common::Vector{Set{Tuple{Int,Int}}},
        neighbours::Vector{Dict{Int,Set{Int}}},
        coms_agents::Vector{Vector{Int}},
        perc_rewire::Float64,
        verbose::Bool=false)

Increases correlation between edges in two selected layers.

**Arguments**
* `layers::Tuple{Int,Int}` IDs of modified layers
* `edges::Vector{Set{Tuple{Int,Int}}}` edges of multilayer network modified in-place to keep state between batches
* `edges_common::Vector{Set{Tuple{Int,Int}}}` edges on common active agents between layers
* `neighbours::Vector{Dict{Int,Set{Int}}}` mapping of neighbours of agents
* `coms_agents::Vector{Vector{Int}}` assignment of agents to communities
* `perc_rewire::Float64` fraction of edges to be used for rewiring attempt
* `verbose::Bool=false` verbose logging flag
"""
function increase_edges_correlation!(layers::Tuple{Int,Int},
    edges::Vector{Set{Tuple{Int,Int}}},
    edges_common::Vector{Set{Tuple{Int,Int}}},
    neighbours::Vector{Dict{Int,Set{Int}}},
    coms_agents::Vector{Vector{Int}},
    perc_rewire::Float64,
    verbose::Bool=false)
    n_rewire = ceil(Int, perc_rewire * minimum(length.(edges_common)))
    adjusted = false
    for _ in 1:n_rewire
        primary, secondary = shuffle([1, 2])
        layer_sec = layers[secondary]
        edges_sec = edges_common[secondary]
        nbs_sec = neighbours[secondary]
        edge = rand(edges_common[primary])
        edge in edges_sec && continue # edge already exists
        u, v = edge
        mult = coms_agents[layer_sec][u] == coms_agents[layer_sec][v] ? 1 : -1
        (!haskey(nbs_sec, mult * u) || length(nbs_sec[mult*u]) == 0) && continue
        (!haskey(nbs_sec, mult * v) || length(nbs_sec[mult*v]) == 0) && continue #no neighbours in community or background
        u_prim = rand(nbs_sec[mult*u])
        v_prim = rand(nbs_sec[mult*v])

        (mult == -1) && (coms_agents[layer_sec][u_prim] == coms_agents[layer_sec][v_prim]) && continue # neighbours from the same comm while u,v from backgorund
        extrema((v_prim, u_prim)) in edges_sec && continue #edge already exists
        length(Set([u, v, u_prim, v_prim])) != 4 && continue #u,v,u_prim,v_prim are 4 distinct nodes
        verbose && println("Wiring $(u_prim)-$(v_prim) and $(u)-$(v). Breaking $(u)-$(u_prim) and $(v)-$(v_prim)")
        adjusted = true
        #Push new edges and neighbours
        @assert !((u, v) in edges[layer_sec])
        @assert !(extrema((u_prim, v_prim)) in edges[layer_sec])
        push!(edges_sec, (u, v), extrema((u_prim, v_prim)))
        push!(edges[layer_sec], (u, v), extrema((u_prim, v_prim)))
        push!(nbs_sec[mult*u], v)
        push!(nbs_sec[mult*v], u)
        push!(nbs_sec[mult*u_prim], v_prim)
        push!(nbs_sec[mult*v_prim], u_prim)
        #Delete edges and neighbours
        @assert extrema((u, u_prim)) in edges[layer_sec]
        @assert extrema((v, v_prim)) in edges[layer_sec]
        delete!(edges_sec, extrema((u, u_prim)))
        delete!(edges_sec, extrema((v, v_prim)))
        delete!(edges[layer_sec], extrema((u, u_prim)))
        delete!(edges[layer_sec], extrema((v, v_prim)))
        delete!(nbs_sec[mult*u], u_prim)
        delete!(nbs_sec[mult*v], v_prim)
        delete!(nbs_sec[mult*u_prim], u)
        delete!(nbs_sec[mult*v_prim], v)
    end
    return adjusted
end

"""
    decrease_edges_correlation!(layers::Tuple{Int,Int},
        edges::Vector{Set{Tuple{Int,Int}}},
        edges_common::Vector{Set{Tuple{Int,Int}}},
        neighbours::Vector{Dict{Int,Set{Int}}},
        coms_agents::Vector{Vector{Int}},
        perc_rewire::Float64,
        verbose::Bool=false)

Decreases correlation between edges in two selected layers.

**Arguments**
* `layers::Tuple{Int,Int}` IDs of modified layers
* `edges::Vector{Set{Tuple{Int,Int}}}` edges of multilayer network modified in-place to keep state between batches
* `edges_common::Vector{Set{Tuple{Int,Int}}}` edges on common active agents between layers
* `edges_in_coms::Vector{Dict{Int,Set{Tuple{Int,Int}}}}` mapping of edges withing communities and background graph
* `coms_agents::Vector{Vector{Int}}` assignment of agents to communities
* `perc_rewire::Float64` fraction of edges to be used for rewiring attempt
* `verbose::Bool=false` verbose logging flag
"""
function decrease_edges_correlation!(layers::Tuple{Int,Int},
    edges::Vector{Set{Tuple{Int,Int}}},
    edges_common::Vector{Set{Tuple{Int,Int}}},
    edges_in_coms::Vector{Dict{Int,Set{Tuple{Int,Int}}}},
    coms_agents::Vector{Vector{Int}},
    perc_rewire::Float64,
    verbose::Bool=false)
    n_rewire = ceil(Int, perc_rewire * minimum(length.(edges_common)))
    adjusted = false
    for _ in 1:n_rewire
        primary, secondary = shuffle([1, 2])
        layer_sec = layers[secondary]
        edges_sec = edges_common[secondary]
        coms_edges_sec = edges_in_coms[secondary]
        common_edges = intersect(edges_common[primary], edges_sec)
        length(common_edges) == 0 && break
        edge = rand(common_edges)
        u, v = edge
        u_com = coms_agents[layer_sec][u]
        v_com = coms_agents[layer_sec][v]
        com = u_com == v_com ? u_com : -1
        u_prim, v_prim = rand(coms_edges_sec[com])
        length(Set([u, v, u_prim, v_prim])) != 4 && continue #u,v,u_prim,v_prim are 4 distinct nodes
        extrema((u, u_prim)) in edges_sec && continue #edge already exists
        extrema((v, v_prim)) in edges_sec && continue #edge already exists
        verbose && println("Wiring $(u)-$(u_prim) and $(v)-$(v_prim). Breaking $(u)-$(v) and $(u_prim)-$(v_prim)")
        adjusted = true
        #Push new edges and neighbours
        new_edge1 = extrema((u, u_prim))
        new_edge2 = extrema((v, v_prim))
        push!(edges_sec, new_edge1, new_edge2)
        push!(edges[layers[secondary]], new_edge1, new_edge2)
        push!(coms_edges_sec[com], new_edge1, new_edge2)
        #Delete edges and neighbours
        delete!(edges_sec, (u, v))
        delete!(edges_sec, (u_prim, v_prim))
        delete!(edges[layers[secondary]], (u, v))
        delete!(edges[layers[secondary]], (u_prim, v_prim))
        delete!(coms_edges_sec[com], (u, v))
        delete!(coms_edges_sec[com], (u_prim, v_prim))
    end
    return adjusted
end

"""
    adjust_edges_correlation(cfg::MLNConfig,
        edges::Vector{Set{Tuple{Int,Int}}},
        coms::Vector{Vector{Int}},
        active_nodes::Vector{Vector{Int}},
        verbose::Bool=false,
        save_cor_change_to_file::Bool=false)

Run `cfg.t` rewiring batches to match edges correlation between layers with desired `cfg.edges_cor_matrix` correlation matrix.

**Arguments**
* `cfg::MLNConfig` MLNABCD configuration
* `edges::Vector{Set{Tuple{Int,Int}}}` agent-based edges of multilayer network before rewiring
* `coms::Vector{Vector{Int}}` assignment of nodes to communities in each layer
* `active_nodes::Vector{Vector{Int}}` IDs of active agents
* `verbose::Bool=false` verbose logging flag
* `save_cor_change_to_file::Bool=false` correlation changes file logging flag
"""
function adjust_edges_correlation(cfg::MLNConfig,
    edges::Vector{Set{Tuple{Int,Int}}},
    coms::Vector{Vector{Int}},
    active_nodes::Vector{Vector{Int}},
    verbose::Bool=false,
    save_cor_change_to_file::Bool=false)

    cfg.skip_edges_correlation && return edges
    best_edges = working_edges = deepcopy(edges)
    common_agents = common_agents_dict(active_nodes)
    edges_common_agents = common_agents_edges(working_edges, common_agents)
    edges_cor, cor_diff, cor_dist = calculate_edges_cor(edges_common_agents, cfg.edges_cor_matrix)
    if save_cor_change_to_file
        now_str = replace(string(now()), r"-|\.|:|T" => "_")
        io = open("edges_correlation_$(now_str).log", "a")
        line = ""
        for i in 1:cfg.l
            for j in (i+1):cfg.l
                line *= string(i) * "-" * string(j) * ","
            end
        end
        println(io, line * "l2_distance")
        write_cor_diff_to_file(io, cor_diff, cor_dist)
    end
    best_cor_dist = Inf
    for b in 1:cfg.t
        layers, direction = pick_cor_to_improve(cor_diff, "weighted")
        verbose && println("Adjusting correlation for layers:", layers, " Direction:", direction)
        i, j = layers
        edges_common = edges_common_agents[(i, j)]
        if direction == "increase"
            neighbours = map_neighbours.(edges_common, coms[[i, j]])
            adjusted = increase_edges_correlation!(layers, working_edges, edges_common, neighbours, coms, cfg.eps, verbose)
        elseif direction == "decrease"
            edges_in_communities = map_communities.(edges_common, coms[[i, j]])
            adjusted = decrease_edges_correlation!(layers, working_edges, edges_common, edges_in_communities, coms, cfg.eps, verbose)
        else
            continue
        end
        if adjusted
            edges_common_agents = common_agents_edges(working_edges, common_agents)
            edges_cor, cor_diff, cor_dist = calculate_edges_cor(edges_common_agents, cfg.edges_cor_matrix)
        end
        if b > 100 && cor_dist < best_cor_dist
            best_cor_dist = cor_dist
            best_edges = deepcopy(working_edges)
            verbose && println("Saved new best edges at batch:", b, " with best distance:", best_cor_dist)
        end
        save_cor_change_to_file && write_cor_diff_to_file(io, cor_diff, cor_dist)
    end
    save_cor_change_to_file && close(io)
    return best_edges
end