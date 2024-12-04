
function sample_ranking(vs, sigma)
    n_vs = length(vs)
    ranks = sortperm([randn() * sigma + v / n_vs for v in vs]) # N(v/n,sigma) for each v
    return ranks
end

function find_ranking(vs, rho, verbose=false, limit=20)
    best_diff = Inf
    best_ranking = []
    #special cases
    if abs(rho) == 1
        best_ranking = collect(vs)
    elseif rho == 0
        for _ in 1:limit
            tmp_ranking = shuffle(vs)
            tmp_sampled_rho = corspearman(vs, tmp_ranking)
            if abs(tmp_sampled_rho - rho) < best_diff
                best_ranking = tmp_ranking
                best_diff = abs(tmp_sampled_rho - rho)
            end
        end
    else
        rho_pos = abs(rho)
        sigma_rho = readdlm("../src/sigma_rho.csv", ',')
        idx = argmin(abs.(sigma_rho[:, 2] .- rho_pos))
        sigma = sigma_rho[idx, 1]
        for _ in 1:limit
            tmp_ranking = sample_ranking(vs, sigma)
            tmp_sampled_rho = corspearman(vs, tmp_ranking)
            if abs(tmp_sampled_rho - rho_pos) < best_diff
                best_ranking = tmp_ranking
                best_diff = abs(tmp_sampled_rho - rho_pos)
            end
        end
    end
    if rho < 0
        best_ranking = reverse(best_ranking)
    end
    best_r = corspearman(vs, best_ranking)
    diff = abs(best_r - rho)
    verbose && println("Input ρ=$(rho) | Generated ρ=$(best_r) | Difference=$(diff)")
    return (rho=best_r, ranking=best_ranking, diff=diff)
end

function degrees_correlation(
    n::Int,
    degs::Vector{Vector{Int}},
    rhos::Vector{Float64},
    active_nodes::Vector{Vector{Int}},
    verbose::Bool=false)
    @assert length(degs) == length(rhos) "Length of correlation vector must be the same as generated degree sequences"
    @assert length.(degs) == length.(active_nodes) "Size of degree sequence must match number of active nodes in each layer"
    vs = 1:n
    orderings = [find_ranking(vs, rho, verbose).ranking for rho in rhos]
    degs_correlated = Vector{Vector{Int}}()
    for i in eachindex(degs)
        degs_layer = zeros(Int, length(degs[i]))
        degs_sorted = sort(degs[i], rev=true)
        if length(active_nodes[i]) == n
            adjusted_ordering = orderings[i]
        else
            mapping = Dict(zip(active_nodes[i],1:length(active_nodes[i])))
            adjusted_ordering = [mapping[e] for e in orderings[i] if e in Set(active_nodes[i])]
        end
        for j in eachindex(degs_sorted)
            degs_layer[adjusted_ordering[j]] = degs_sorted[j]
        end
        push!(degs_correlated, copy(degs_layer))
    end
    return degs_correlated
end

struct MLNABCDParams
    w::Vector{Int}
    s::Vector{Int}
    ξ::Float64
    μ::Nothing
    isCL::Bool
    islocal::Bool
    hasoutliers::Bool

    function MLNABCDParams(w, s, ξ)
        length(w) == sum(s) || throw(ArgumentError("inconsistent data"))
        0 ≤ ξ ≤ 1 || throw(ArgumentError("inconsistent data ξ"))

        new(w, s, ξ, nothing, false, false, false)
    end
end