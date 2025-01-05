
function sample_points(n::Int, d::Int)
    points = randn(n, d)
    points ./= sqrt.(sum(x -> x^2, points, dims=2))
    points .*= rand(n) .^ (1 / d)
    return points
end

function assign_points(x::Matrix{Float64}, c::Vector{Int})
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

function shuffle_communities(r::Float64, a::Vector{Vector{Int}})
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

function communities_correlation(n::Int,
    coms::Vector{Vector{Int}},
    rs::Vector{Float64},
    active_nodes::Vector{Vector{Int}})

    x = sample_points(n, 2)
    coms_correlated = Vector{Vector{Int}}()
    for i in eachindex(coms)
        x_active = x[active_nodes[i], :]
        a = assign_points(x_active, coms[i])
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

function generate_communities(cfg::MLNConfig, active_nodes::Vector{Vector{Int}})
    coms_sizes = ABCDGraphGenerator.sample_communities.(cfg.betas, cfg.c_mins, cfg.c_maxs, cfg.ns, Ref(cfg.c_max_iter))
    return coms_sizes, communities_correlation(cfg.n, coms_sizes, cfg.rs, active_nodes)
end