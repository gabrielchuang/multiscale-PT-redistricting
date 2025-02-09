import Pkg
#codepath = joinpath("/home/postdoc/gtc17/MSMS_copy", "MultiScaleMapSampler") 
codepath = joinpath("../../") 
Pkg.activate("mergeSplit"; shared=true)
Pkg.develop(path=codepath)

using Revise
using MultiScaleMapSampler
using JSON, CSV, Test, Random, RandomNumbers, TickTock, Dates

## ==== LOAD THE GRAPH ==== 
n = 16
G_json = joinpath("..","test_graphs", "$(n)x$(n)pct.json")
G_node_data = Set(["county", "pct", "pop", "area", "border_length", "id"])
G_base_graph = BaseGraph(G_json, "pop", inc_node_data=G_node_data,
	area_col="area",node_border_col="border_length", edge_perimeter_col="length")
for k in G_base_graph.node_attributes 
	k["county_and_pct"] = string(k["county"], ",",k["pct"])
end 

G = MultiLevelGraph(G_base_graph, ["county_and_pct"]);
oriented_nbrs = JSON.parsefile("../test_graphs/$(n)x$(n)pct_oriented_nbrs.json")
oriented_nbrs = Dict{GlobalID, Vector{GlobalID}}([parse(Int, k) => v for (k,v) in oriented_nbrs]);

## ==== LOAD THE HIERARCHY ====
filename_prefix = "../test_graphs/$(n)x$(n)TH_5/tempered_hierarchy_"
TH = TemperingHierarchy(G, filename_prefix, 9);
num_levels = length(TH.graphs_by_level)
num_dists = 4

rngseed = 15300
rng = PCG.PCGStateOneseq(UInt64, rngseed)
suffix = "9L_F"

SNF_step_sizes = [30 for _=1:num_levels]

total_pop = n*n; ideal_district_pop = total_pop/num_dists
max_node_sizes = [maximum(filter(y->y != n*n, [x["pop"] for x in values(TH.graphs_by_level[i].node_attributes)])) for i =1:num_levels]
num_nodes_by_level = [length(finest_level_nodes(x)) for x in TH.graphs_by_level] 
pop_devs = (max_node_sizes .+ 1.5) ./ ideal_district_pop
@show max_node_sizes, pop_devs
ideal_pop = get_node_population(G, G.root.global_ID) / num_dists
pop_mins_maxs = [(max(0.01, (ideal_pop*(1-pop_dev))), ideal_pop*(1+pop_dev)) for pop_dev in pop_devs]


## ==== SET THE MEASURE AT EACH LEVEL OF THE HIERARCHY ==== 
iso_weights = [0.3 for _=1:num_levels]
SF_weights = [-1.0 for _=1:num_levels]
betas = [0.1 for _=1:num_levels]

constraints = []
measures = []
temper_measures = []
for j=1:num_levels 
	constraint = initialize_constraints() 
	add_constraint!(constraint, PopulationConstraint(pop_mins_maxs[j][1], pop_mins_maxs[j][2]))
	push!(constraints, constraint)

	measure = Measure(0.0)
	push_energy!(measure, compactnessEnergy(nothing), iso_weights[j], "compactness")
	push!(measures, measure)
	push_energy!(measure, spanningForestEnergy(0), iso_weights[j], "spanning forest count")
	push!(measures, measure)

	temper_measure = Measure(0.0)
	push!(temper_measures, temper_measure)
end   

rank = 1
processor_to_part_ID = [collect(1:num_levels)]

# ==== OUTPUT SETUP  ====
writers::Vector{Union{SimpleWriter, Nothing}} = [nothing for _=1:num_levels]
output_file_path_1 = "../output/$(n)x$(n)_districts/"*suffix*".jsonl"
writers[1] = SimpleWriter(measures[1], constraints[1], num_dists, output_file_path_1)

part_ID_to_temp::Dict{Int, Int} = Dict([x => x for x in 1:num_levels])
my_partitions = Dict{Int, MultiLevelPartition}() #part_ID to partition 

# ==== INITIALIZE PARTITIONS AT EACH TEMPERATURE ====
for part_ID in 1:num_levels
	@show part_ID
	t = part_ID_to_temp[part_ID]
	partition = MultiLevelPartition(TH.graphs_by_level[t], constraints[t], num_dists; rng=rng, max_attempts=100)
	populate_extras!(TH.graphs_by_level[t])
	partition.extensions["iso"] = []; partition.extensions["aps"] = []; partition.extensions["pop"] = []
	partition.extensions["temperature"] = []
	partition.extensions["partition_ID"] = part_ID
	my_partitions[part_ID] = partition

	run_metropolis_hastings_SNF!(my_partitions[part_ID], measures[t], temper_measures[t], 
		100, constraints[t], rng, writer=writers[t], output_freq=100)
end

inner_steps = 10000 # NUMBER OF SWAP ATTEMPTS TO MAKE 

println("total SNF steps: ", outer_steps*inner_steps*SNF_step_sizes[end], "; total PT flip attempts: ", outer_steps*inner_steps)

total_steps_taken = Dict([x => 0 for x in 1:num_levels])

# ==== LOGGING ==== 
no_merges = [0 for i=1:num_levels]
no_splits = [0 for i=1:num_levels]
almosts = [0 for i=1:num_levels] 
succs = [0 for i=1:num_levels]

# ==== RUN THE SAMPLER ==== 

for i=1:inner_steps 
	if mod(i, 500) == 0 
		@show i, now()
		flush(stdout)
	end

	for part_ID in 1:num_levels
		t = part_ID_to_temp[part_ID] 
		SNF_step_size = SNF_step_sizes[t]
		run_metropolis_hastings_SNF!(my_partitions[part_ID], measures[t], temper_measures[t], 
			(total_steps_taken[part_ID], total_steps_taken[part_ID]+SNF_step_size), constraints[t], rng, writer=writers[t], output_freq=25)
		total_steps_taken[part_ID] += SNF_step_size
		push!(my_partitions[part_ID].extensions["temperature"], t)
	end 
	PT_swap!(TH, my_partitions, part_ID_to_temp, num_levels, pop_mins_maxs, measures, i, rng, betas, 
		(no_merges, no_splits, almosts, succs))
end

@show no_merges,no_splits,almosts,succs