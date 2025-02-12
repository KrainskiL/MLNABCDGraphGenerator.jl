
"""
    sample_points(n::Int, d::Int)::Matrix{Float64}

Samples normalized points in `d`-dimensional space.

**Arguments**
* `n::Int` number of points
* `d::Int` number of dimensions
"""
function sample_points(n::Int, d::Int)::Matrix{Float64}
    points = randn(n, d)
    points ./= sqrt.(sum(x -> x^2, points, dims=2))
    points .*= rand(n) .^ (1 / d)
    return points
end

"""
    assign_points(x::Matrix{Float64}, c::Vector{Int})::Vector{Vector{Int}}

Assigns points to communities based on the proximity in Euclidian space.

**Arguments**
* `x::Matrix{Float64}` coordinates of points
* `c::Vector{Int}` sizes of communities
"""
function assign_points(x::Matrix{Float64}, c::Vector{Int})::Vector{Vector{Int}}
    @assert ndims(x) == 2
    @assert sum(c) == size(x, 1)
    c = shuffle(c)
    x = copy(x)
    all_idxs = collect(1:size(x, 1))
    dist = vec(sum(x -> x^2, x, dims=2))
    res = Vector{Int}[]
    for com in c
        ind = argmax(dist)
        ref = x[ind, :]
        dist_c = [sum(x -> x^2, (r - ref)) for r in eachrow(x)]
        idxs = partialsortperm(dist_c, 1:com)
        push!(res, all_idxs[idxs])
        to_keep = setdiff(1:size(x, 1), idxs)
        x = x[to_keep, :]
        dist = dist[to_keep]
        all_idxs = all_idxs[to_keep]
    end
    @assert size(x, 1) == 0
    @assert length(all_idxs) == 0
    @assert length(dist) == 0
    @assert sort(union(res...)) == 1:sum(c)
    return res
end

"""
    shuffle_communities(r::Float64, a::Vector{Vector{Int}})::Vector{Vector{Int}}

Shuffles assignment of nodes to communities to introduce correlation `r` between the layer and reference layer

**Arguments**
* `r::Float64` desired correlation between the layer and reference layer
* `a::Vector{Vector{Int}}` assignment of nodes to communities
"""
function shuffle_communities(r::Float64, a::Vector{Vector{Int}})::Vector{Vector{Int}}
    shuffle_nodes = []
    for i in 1:length(a)
        for j in 1:length(a[i])
            if rand() < 1 - r
                push!(shuffle_nodes, (i, j))
            end
        end
    end

    shuffled_a = deepcopy(a)
    if length(shuffle_nodes) > 1
        shuffle!(shuffle_nodes)
        for i in 1:length(shuffle_nodes)-1
            com1, pos1 = shuffle_nodes[i]
            com2, pos2 = shuffle_nodes[i+1]
            tmp_node = shuffled_a[com1][pos1]
            shuffled_a[com1][pos1] = shuffled_a[com2][pos2]
            shuffled_a[com2][pos2] = tmp_node
        end
    end
    return shuffled_a
end

"""
    communities_correlation(
        n::Int,
        coms_sizes::Vector{Vector{Int}},
        rs::Vector{Float64},
        active_nodes::Vector{Vector{Int}},
        d::Int)
        ::Vector{Vector{Int}}

Shuffles communities assignment to match correlations `rs` between network layers and hidden reference layer.

**Arguments**
* `n::Int` number of agents
* `coms_sizes::Vector{Vector{Int}}` sizes of communities
* `rs::Vector{Float64}` vector of desired correlations between each network layer and reference layer
* `active_nodes::Vector{Vector{Int}}` IDs of active nodes
* `d::Int` number of dimensions in reference space
"""
function communities_correlation(
    n::Int,
    coms_sizes::Vector{Vector{Int}},
    rs::Vector{Float64},
    active_nodes::Vector{Vector{Int}},
    d::Int)::Vector{Vector{Int}}

    x = sample_points(n, d)
    coms_correlated = Vector{Vector{Int}}()
    for i in eachindex(coms_sizes)
        x_active = x[active_nodes[i], :]
        a = assign_points(x_active, coms_sizes[i])
        shuffled_a = shuffle_communities(rs[i], a)
        flattened = zeros(Int, length(active_nodes[i]))
        for (j, com) in enumerate(shuffled_a)
            for e in com
                flattened[e] = j
            end
        end
        push!(coms_correlated, flattened)
    end
    return coms_correlated
end

"""
    generate_communities(cfg::MLNConfig, active_nodes::Vector{Vector{Int}})

Samples sequences of community sizes from truncated power-law distribution.
Then, assigns nodes to communities and shuffles the assignment to correlate layers.
Returns sizes of communities and assignment of nodes to communities.

**Arguments**
* `cfg::MLNConfig` MLNABCD configuration
* `active_nodes::Vector{Vector{Int}}` IDs of active nodes
"""
function generate_communities(cfg::MLNConfig, active_nodes::Vector{Vector{Int}})
    coms_sizes = ABCDGraphGenerator.sample_communities.(cfg.betas, cfg.c_mins, cfg.c_maxs, cfg.ns, Ref(cfg.c_max_iter))
    return coms_sizes, communities_correlation(cfg.n, coms_sizes, cfg.rs, active_nodes, cfg.d)
end