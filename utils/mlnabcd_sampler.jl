using MLNABCDGraphGenerator

@info "Usage: julia mlnabcd_sampler.jl config_filename.toml"
@info "For the syntax of config_filename see example_config.toml file"

filename = ARGS[1]
config = MLNABCDGraphGenerator.parse_config(filename)
#Active nodes
active_nodes = MLNABCDGraphGenerator.generate_active_nodes(config)
#Degree Sequences
degrees = MLNABCDGraphGenerator.generate_degrees(config, active_nodes, false)
#Sizes of communities
com_sizes, coms = MLNABCDGraphGenerator.generate_communities(config, active_nodes)
#Generate ABCD graphs
edges = MLNABCDGraphGenerator.generate_abcd(config, degrees, com_sizes, coms)
#Map nodes and communities to agents
edges = MLNABCDGraphGenerator.map_edges_to_agents(edges, active_nodes)
coms = MLNABCDGraphGenerator.map_communities_to_agents(config.n, coms, active_nodes)
#Adjust edges correlation
edges_rewired = MLNABCDGraphGenerator.adjust_edges_correlation(config, edges, coms, active_nodes, false, false)
#Save edges to file
MLNABCDGraphGenerator.write_edges(config, edges_rewired)
#Save communities to file
MLNABCDGraphGenerator.write_edges(config, coms)
