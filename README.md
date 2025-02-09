# Multiscale Parallel Tempering Sampler for Redistricting

This codebase implements the sampling algorithm from [Multiscale Parallel Tempering for Fast Sampling on Redistricting Plans](https://arxiv.org/abs/2401.17455). 

## Files
The implementation is found in `src/`. Examples of scripts, as well as sample graphs for square grids (4x4, 8x8, 16x16) and Connecticut, can be found in `examples/`. The interested reader should follow the following sequence: 

1. Sample graphs are found in `examples/test_graphs/`.
2. To generate a hierarchy on a base graph, see the scripts in `examples/hierarchy_creation`. Sample hierarchies (such as `16x16TH_6`) can be found in `examples/test_graphs`. 
3. A list of oriented neighbors - that is, for each node, a clockwise ordered list of neighbors - is required, and can be generated using `examples/preprocessing/generate_oriented_nbrs.ipynb`. 
4. Given the base graph, tempering hierarchy, and oriented neighbors, the algorithm is run using the scripts in `examples/runscripts`, generating files in some `output` directory. 
	- For single-threaded running, see `CT_hierarchy_swapping.jl` or `small_grid_hierarchy_swapping.jl`.  
	- For MPI multi-processor usage, see `CT_hierarchy_swapping_multiproc.jl`. Note that this file is substantially longer, because of MPI calls interleaved in the code. 
5. One can compute observables (e.g., isoperimetric scores, spanning forest counts, vote share marginals) on the output districts using the scripts in `examples/post_processing`. 

The Connecticut sample data files are based on original data from the [Redistricting Data Hub](https://redistrictingdatahub.org/). 