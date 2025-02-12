
"""
    sample_ranking(vs::Vector{Int}, n::Int, sigma::AbstractFloat)::Vector{Vector{Int}}

Returns ranking of samples of random variables N(v/n,sigma).

**Arguments**
* `vs::Vector{Int}` vector with IDs of active agents
* `n::Int` number of agents
* `sigma::AbstractFloat` standard deviation of normal distribution N(v/n,sigma)
"""
function sample_ranking(vs::Vector{Int}, n::Int, sigma::AbstractFloat)::Vector{Int}
    ranks = sortperm([randn() * sigma + v / n for v in vs])
    return ranks
end

"""
    find_ranking(ns::Int, n::Int, tau::AbstractFloat, verbose::Bool=false, limit::Int=20)

Finds ordering of labels of agents which has `tau` Kendall-tau correlation with ordered labels sequence.

**Arguments**
* `ns::Int` number of active agents
* `n::Int` number of agents
* `tau::AbstractFloat` desired Kendall-tau correlation value
* `verbose::Bool=false` verbose logging flag
* `limit::Int=20` number of rankings sampled to find the best correlation match
"""
function find_ranking(ns::Int, n::Int, tau::AbstractFloat, verbose::Bool=false, limit::Int=20)
    vs = collect(1:ns)
    best_diff = Inf
    best_ranking = []
    #special cases
    if abs(tau) == 1
        best_ranking = vs
    elseif tau == 0
        for _ in 1:limit
            tmp_ranking = shuffle(vs)
            tmp_sampled_tau = corkendall(vs, tmp_ranking)
            if abs(tmp_sampled_tau - tau) < best_diff
                best_ranking = tmp_ranking
                best_diff = abs(tmp_sampled_tau - tau)
            end
        end
    else
        tau_pos = abs(tau)
        sigma_tau = readdlm(joinpath(@__DIR__, "sigma_tau.csv"), ',')
        idx = argmin(abs.(sigma_tau[:, 2] .- tau_pos))
        sigma = sigma_tau[idx, 1]
        for _ in 1:limit
            tmp_ranking = sample_ranking(vs, n, sigma)
            tmp_sampled_tau = corkendall(vs, tmp_ranking)
            if abs(tmp_sampled_tau - tau_pos) < best_diff
                best_ranking = tmp_ranking
                best_diff = abs(tmp_sampled_tau - tau_pos)
            end
        end
    end
    if tau < 0
        best_ranking = reverse(best_ranking)
    end
    best_r = corkendall(vs, best_ranking)
    diff = abs(best_r - tau)
    verbose && println("Input τ=$(tau) | Generated ρ=$(best_r) | Difference=$(diff)")
    return (tau=best_r, ranking=best_ranking, diff=diff)
end

"""
    degrees_correlation(
        cfg::MLNConfig,
        degs::Vector{Vector{Int}},
        active_nodes::Vector{Vector{Int}},
        verbose::Bool=false)
        ::Vector{Vector{Int}}

Creates correlated sequences of degrees for all layers of the network.

**Arguments**
* `cfg::MLNConfig` MLNABCD configuration
* `degs::Vector{Vector{Int}}` sorted (input) sequences of degrees
* `active_nodes::Vector{Vector{Int}}` IDs of active nodes
* `verbose::Bool=false` verbose logging flag
"""
function degrees_correlation(
    cfg::MLNConfig,
    degs::Vector{Vector{Int}},
    active_nodes::Vector{Vector{Int}},
    verbose::Bool=false)::Vector{Vector{Int}}
    ns = length.(active_nodes)
    @assert length(degs) == length(cfg.taus) "Length of correlation vector must be the same as generated degree sequences"
    @assert length.(degs) == ns "Size of degree sequence must match number of active nodes in each layer"
    orderings = find_ranking.(ns, Ref(cfg.n), cfg.taus, verbose)
    orderings = getfield.(orderings, :ranking)
    degs_correlated = deepcopy(degs)
    for i in eachindex(degs)
        degs_correlated[i][orderings[i]] .= degs[i]
    end
    return degs_correlated
end

"""
    generate_degrees(cfg::MLNConfig,
        active_nodes::Vector{Vector{Int}},
        verbose::Bool=false)
        ::Vector{Vector{Int}}

Samples sequences of degrees from truncated power-law distribution and reorder the sequences to correlate layers.

**Arguments**
* `cfg::MLNConfig` MLNABCD configuration
* `active_nodes::Vector{Vector{Int}}` IDs of active nodes
* `verbose::Bool=false` verbose logging flag
"""
function generate_degrees(
    cfg::MLNConfig,
    active_nodes::Vector{Vector{Int}},
    verbose::Bool=false)::Vector{Vector{Int}}
    degs = ABCDGraphGenerator.sample_degrees.(cfg.gammas, cfg.d_mins, cfg.d_maxs, cfg.ns, Ref(cfg.d_max_iter))
    return degrees_correlation(cfg, degs, active_nodes, verbose)
end