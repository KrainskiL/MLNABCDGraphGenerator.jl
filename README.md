# MLNABCDGraphGenerator.jl
## Installation
Run the following commands in your Julia REPL:
```julia
using Pkg
Pkg.add("https://github.com/KrainskiL/MLNABCDGraphGenerator.jl")

```

## Usage
The main file intended to be used is `mlnabcd_sampler.jl`.
The script requires a configuration file with global parameters (see `utils/example_config.toml`).
Example contents of the file:
```
seed = "42"                         # RNG seed, use "" for no seeding
n = "1000"                          # number of agents
R = "0.2"                           # desired edges correlation
layer_params = "layer_params.csv"   # file with parameters for each layer
d_max_iter = "1000"                 # maximum number of iterations for sampling degrees
c_max_iter = "1000"                 # maximum number of iterations for sampling cluster sizes
t = "1000"                          # number of batches for edge rewiring
e = "0.1"                           # percent of edges to be rewired in each rewiring batch
```
Layer-specific parameters are provided in file referenced under layer_params key.
Example layers specification formatted as Markdown table (see also `utils/layer_params.csv`):

| q    | rho  | r    | gamma | delta | Delta | beta | s | S  | xi  |
|------|------|------|-------|-------|-------|------|---|----|-----|
| 1    | 1    | 1    | 2,5   | 2     | 25    | 1,5  | 8 | 32 | 0,2 |
| 0,75 | 0,75 | 0,75 | 2,5   | 2     | 25    | 1,5  | 8 | 32 | 0,2 |
| 0,5  | 0,5  | 0,5  | 2,5   | 2     | 20    | 1,7  | 8 | 32 | 0,2 |
| 0,25 | 0,25 | 0,25 | 2,5   | 2     | 20    | 1,7  | 8 | 32 | 0,1 |

The two files contain all parameters required to generate an MLNABCD graph.

Here is an output from an example session using CLI:
```
$ julia mlnabcd_sampler.jl example_config.toml
[ Info: Usage: julia mlnabcd_sampler.jl config_filename
[ Info: For the syntax of config_filename see example_config.toml file
```
After the program terminates two files, `communities.dat` and `edges.dat` are created in the working directory.

Structure of the `communities.dat` is as follows:
```
community_number layer_number
```
Each layer contains assignments for all `n` agents. Inactive agents have are assigned to community 0.

Structure of the `edges.dat` is as follows:
```
agent1_label agent2_label layer_number
```
Edges are based on labels of active agents.
