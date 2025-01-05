
function map_edges_to_agents(edges::Vector{Set{Tuple{Int64,Int64}}}, active_nodes::Vector{Vector{Int64}})
    edges_agents = Vector{Set{Tuple{Int64,Int64}}}()
    for (i, layer) in enumerate(edges)
        push!(edges_agents, Set([(active_nodes[i][edge[1]], active_nodes[i][edge[2]]) for edge in layer]))
    end
    return edges_agents
end

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

function common_agents_edges(edges_agents::Vector{Set{Tuple{Int64,Int64}}}, active_nodes::Vector{Vector{Int64}})
    edges_on_common_agents_dict = Dict()
    for i in eachindex(edges_agents)
        for j in (i+1):length(edges_agents)
            common_agents = intersect(Set(active_nodes[i]), Set(active_nodes[j]))
            edges_i = Set([e for e in edges_agents[i] if (e[1] in common_agents) && (e[2] in common_agents)])
            edges_j = Set([e for e in edges_agents[j] if (e[1] in common_agents) && (e[2] in common_agents)])
            edges_on_common_agents_dict[(i, j)] = (copy(edges_i), copy(edges_j))
        end
    end
    return edges_on_common_agents_dict
end

function map_neighbours(edges, coms)
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

function map_communities(edges, coms)
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
function neighbours_common_agents_edges(edges_common_agents, coms_agents::Vector{Vector{Int64}})
    nbhs_edges_on_common_agents_dict = Dict()
    for (k, v) in edges_common_agents
        nbhs_edges_on_common_agents_dict[k] = (map_neighbours(v[1], coms_agents[k[1]]), map_neighbours(v[2], coms_agents[k[2]]))
    end
    return nbhs_edges_on_common_agents_dict
end

function edges_in_communities_common_agents_edges(edges_common_agents, coms_agents::Vector{Vector{Int64}})
    edges_in_coms = Dict()
    for (k, v) in edges_common_agents
        edges_in_coms[k] = (map_communities(v[1], coms_agents[k[1]]), map_communities(v[2], coms_agents[k[2]]))
    end
    return edges_in_coms
end

function calculate_edges_cor(edges_common_agents, rounded=false)
    current_rs = []
    for (k, v) in edges_common_agents
        e1, e2 = v
        r = length(intersect(e1, e2)) / (min(length(e1), length(e2)))
        if rounded
            r = round(r, digits=4)
        end
        push!(current_rs, (k, r))
    end
    return current_rs
end

function pick_cor_to_improve(cor_vec, desired_cor::Matrix{Float64}, method::String="weighted")
    @assert method in ["weighted", "uniform", "e-greedy"] "Method must bed weighted, uniform or e-greedy, got $(method)"
    cor_diff = []
    for (layers, current_cor) in cor_vec
        diff = current_cor - desired_cor[layers...]
        diff = isnan(diff) ? 0.0 : diff
        push!(cor_diff, (layers, diff))
    end
    cor_diff_values = abs.(getindex.(cor_diff, 2))
    if method == "uniform"
        pick = rand(cor_diff)
    elseif method == "weighted"
        pick = sample(cor_diff, Weights(cor_diff_values))
    else
        ϵ = 0.5
        max_diff_idx = argmax(cor_diff_values)
        other_weight = (1 - ϵ) / (length(cor_diff) - 1)
        weights = fill(other_weight, length(cor_diff))
        weights[max_diff_idx] = ϵ
        pick = sample(cor_diff, Weights(weights))
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

function write_cor_diff_to_file(io, cor_vec, desired_cor::Matrix{Float64})
    cor_diff = [(join(e[1], '-'), e[2] - desired_cor[e[1][1], e[1][2]]) for e in cor_vec]
    println(io, join(getindex.(sort(cor_diff, by=x -> x[1]), 2), ','))
end

function increase_edges_correlation!(layers::Tuple{Int,Int}, edges, edges_common, neighbours, coms_agents, perc_rewire::Float64=0.4, verbose::Bool=false)
    n_rewire = ceil(Int, perc_rewire * min(length(edges[1]), length(edges[2])))
    adjusted = 0
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
        adjusted += 1
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
    adjust_ratio = round(adjusted * 100 / n_rewire, digits=2)
    verbose && println("Adjusted edges percentage:", adjust_ratio)
    return adjust_ratio
end

function decrease_edges_correlation!(layers::Tuple{Int,Int}, edges, edges_common, edges_in_coms, coms_agents, perc_rewire::Float64=0.4, verbose::Bool=false)
    n_rewire = ceil(Int, perc_rewire * min(length(edges[1]), length(edges[2])))
    adjusted = 0
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
        adjusted += 1
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
    adjust_ratio = round(adjusted * 100 / n_rewire, digits=2)
    verbose && println("Adjusted edges percentage:", adjust_ratio)
    return adjust_ratio
end

function adjust_edges_correlation!(cfg::MLNConfig, edges, coms::Vector{Vector{Int}},
    active_nodes::Vector{Vector{Int}},
    verbose::Bool=false, save_cor_change_to_file::Bool=false)

    cfg.skip_edges_correlation && return nothing
    edges_common_agents = common_agents_edges(edges, active_nodes)
    edges_cor = calculate_edges_cor(edges_common_agents)
    avg_adjust_ratio = []
    if save_cor_change_to_file
        now_str = replace(string(now()), r"-|\.|:|T" => "_")
        io = open("edges_correlation_$(now_str).log", "a")
        println(io, join(sort(join.(getindex.(edges_cor, 1), '-')), ','))
        write_cor_diff_to_file(io, edges_cor, cfg.edges_cor_matrix)
    end
    for b in 1:cfg.t
        layers, direction = pick_cor_to_improve(edges_cor, cfg.edges_cor_matrix, "weighted")
        verbose && println("Adjusting correlation for layers:", layers, " Direction:", direction)
        i, j = layers
        common_agents = intersect(Set(active_nodes[i]), Set(active_nodes[j]))
        edges1 = Set([e for e in edges[i] if (e[1] in common_agents) && (e[2] in common_agents)])
        edges2 = Set([e for e in edges[j] if (e[1] in common_agents) && (e[2] in common_agents)])
        neighbours = (map_neighbours(edges1, coms[i]), map_neighbours(edges2, coms[j]))
        edges_in_communities = (map_communities(edges1, coms[i]), map_communities(edges2, coms[j]))
        if direction == "increase"
            adj_ratio = increase_edges_correlation!(layers, edges, (edges1, edges2), neighbours, coms, cfg.eps, verbose)
        elseif direction == "decrease"
            adj_ratio = decrease_edges_correlation!(layers, edges, (edges1, edges2), edges_in_communities, coms, cfg.eps, verbose)
        else
            continue
        end
        push!(avg_adjust_ratio, adj_ratio)
        edges_common_agents = common_agents_edges(edges, active_nodes)
        edges_cor = calculate_edges_cor(edges_common_agents)
        save_cor_change_to_file && write_cor_diff_to_file(io, edges_cor, cfg.edges_cor_matrix)
    end
    verbose && print("Average adjust hit rate:", mean(avg_adjust_ratio))
    save_cor_change_to_file && close(io)
    return nothing
end