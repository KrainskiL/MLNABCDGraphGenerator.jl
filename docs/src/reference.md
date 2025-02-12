Reference
=========

```@meta
CurrentModule = MLNABCDGraphGenerator
```

Utilities
----------------------
```@docs
MLNConfig
parse_config
write_edges
write_communities
```

Phase 1 - Active Nodes
----------------------
```@docs
generate_active_nodes
```

Phase 2 - Degrees Correlation
----------------------
```@docs
sample_ranking
find_ranking
degrees_correlation
generate_degrees
```

Phase 3 - Communities Correlation
----------------------
```@docs
sample_points
assign_points
shuffle_communities
communities_correlation
generate_communities
```

Phase 4 and 5 - Edges generation
----------------------
```@docs
generate_abcd
```

Phase 6 - Edges correlation
----------------------
```@docs
map_edges_to_agents
map_communities_to_agents
common_agents_dict
common_agents_edges
map_neighbours
map_communities
neighbours_common_agents_edges
edges_in_communities_common_agents_edges
calculate_edges_cor
pick_cor_to_improve
write_cor_diff_to_file
increase_edges_correlation!
decrease_edges_correlation!
adjust_edges_correlation
```