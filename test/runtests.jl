using Test
using MLNABCDGraphGenerator

cd("../utils")
config_path = "example_config.toml"
cfg = MLNABCDGraphGenerator.parse_config(config_path)

@testset "config" begin
    @test typeof(MLNABCDGraphGenerator.parse_config(config_path)) == MLNABCDGraphGenerator.MLNConfig
    args = [nothing, 10, "", "example_layer_params.csv", 100, 100, 100, 0.01, 2, "edges.dat", "coms.dat"]
    cfg_no_edges_cor = MLNABCDGraphGenerator.MLNConfig(args...)
    @test cfg_no_edges_cor.skip_edges_correlation
    @test sum(cfg_no_edges_cor.edges_cor_matrix) == 0
    args[3] = "0.5"
    cfg_fixed_edges_cor = MLNABCDGraphGenerator.MLNConfig(args...)
    @test !cfg_fixed_edges_cor.skip_edges_correlation
    @test all(cfg_fixed_edges_cor.edges_cor_matrix .== 0.5)
end

@testset "active_nodes" begin
    active_nodes = MLNABCDGraphGenerator.generate_active_nodes(cfg)
    @test typeof(active_nodes) == Vector{Vector{Int}}
    @test issorted(active_nodes[1])
end

active_nodes = MLNABCDGraphGenerator.generate_active_nodes(cfg)
@testset "degrees" begin
    degrees = MLNABCDGraphGenerator.generate_degrees(cfg, active_nodes, false)
    @test typeof(degrees) == Vector{Vector{Int}}
    @test all(length.(active_nodes) == length.(degrees))
end

@testset "communities" begin
    com_sizes, coms = MLNABCDGraphGenerator.generate_communities(cfg, active_nodes)
    @test typeof(com_sizes) == Vector{Vector{Int}}
    @test typeof(coms) == Vector{Vector{Int}}
    @test all(length.(active_nodes) == length.(coms))
    cfg.d = 1
    com_sizes, coms = MLNABCDGraphGenerator.generate_communities(cfg, active_nodes)
    @test typeof(com_sizes) == Vector{Vector{Int}}
    @test typeof(coms) == Vector{Vector{Int}}
    @test all(length.(active_nodes) == length.(coms))
    cfg.d = 3
    com_sizes, coms = MLNABCDGraphGenerator.generate_communities(cfg, active_nodes)
    @test typeof(com_sizes) == Vector{Vector{Int}}
    @test typeof(coms) == Vector{Vector{Int}}
    @test all(length.(active_nodes) == length.(coms))
end

degrees = MLNABCDGraphGenerator.generate_degrees(cfg, active_nodes, false)
com_sizes, coms = MLNABCDGraphGenerator.generate_communities(cfg, active_nodes)
@testset "abcd" begin
    edges = MLNABCDGraphGenerator.generate_abcd(cfg, degrees, com_sizes, coms)
    @test typeof(edges) == Vector{Set{Tuple{Int,Int}}}
end

edges = MLNABCDGraphGenerator.generate_abcd(cfg, degrees, com_sizes, coms)
@testset "agents_remap" begin
    edges_remapped = MLNABCDGraphGenerator.map_edges_to_agents(edges, active_nodes)
    @test typeof(edges_remapped) == Vector{Set{Tuple{Int,Int}}}
    coms_remapped = MLNABCDGraphGenerator.map_communities_to_agents(cfg.n, coms, active_nodes)
    @test typeof(coms_remapped) == Vector{Vector{Int}}
    @test all(length.(coms_remapped) .== cfg.n)
end

edges = MLNABCDGraphGenerator.map_edges_to_agents(edges, active_nodes)
coms = MLNABCDGraphGenerator.map_communities_to_agents(cfg.n, coms, active_nodes)
@testset "edges_rewiring" begin
    edges_rewired = MLNABCDGraphGenerator.adjust_edges_correlation(cfg, edges, coms, active_nodes, false, false)
    @test typeof(edges_rewired) == Vector{Set{Tuple{Int,Int}}}
    cfg.skip_edges_correlation = true
    edges_rewired = MLNABCDGraphGenerator.adjust_edges_correlation(cfg, edges, coms, active_nodes, false, false)
    @test typeof(edges_rewired) == Vector{Set{Tuple{Int,Int}}}
end