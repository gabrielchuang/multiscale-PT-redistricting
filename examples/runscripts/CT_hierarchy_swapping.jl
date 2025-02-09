import Pkg
#codepath = joinpath("/home/postdoc/gtc17/MSMS_copy", "MultiScaleMapSampler") 
codepath = joinpath("../../") 
Pkg.activate("mergeSplit"; shared=true)
Pkg.develop(path=codepath)

using Revise
using MultiScaleMapSampler
using JSON, CSV, Test, Random, RandomNumbers, TickTock, Dates
#using Plots, StatsPlots, GraphRecipes 

CT_json = joinpath("../test_graphs", "CT_pct20_1indexed.json")
CT_node_data = Set(["COUNTY", "NAME", "POP20", "area", "border_length"])
CT_base_graph = BaseGraph(CT_json, "POP20", inc_node_data=CT_node_data,
	area_col="area",node_border_col="border_length", edge_perimeter_col="length")
for k in CT_base_graph.node_attributes 
	k["county_and_pct"] = string(k["COUNTY"], ",",k["NAME"])
end 

Gct = MultiLevelGraph(CT_base_graph, ["county_and_pct"]);
oriented_nbrs = JSON.parsefile("../test_graphs/CT_pct20_oriented_nbrs.json")
oriented_nbrs = Dict{GlobalID, Vector{GlobalID}}([parse(Int, k) => v for (k,v) in oriented_nbrs]);

filename_prefix = "../test_graphs/CTTH_17/tempered_hierarchy_"
TH = TemperingHierarchy(Gct, filename_prefix, 24);
num_levels = length(TH.graphs_by_level)
num_dists = 5

rngseed = length(ARGS) > 1 ? parse(Int64, ARGS[2]) : 15192
rng = PCG.PCGStateOneseq(UInt64, rngseed)
suffix = ARGS[1]

SNF_step_sizes = [30 for _=1:num_levels]


# =========== MEASURES/CONSTRAINTS SETUP	
num_nodes_by_level = [length(finest_level_nodes(x)) for x in TH.graphs_by_level] 
pop_devs = 2.785 ./ (num_nodes_by_level ./ num_dists)

iso_weights = 1 + 0.9 .- sqrt.(1 .+ (collect(1:num_levels) .* 1.616 ./ num_levels) .^ 2) #hyperbola from 0.9 
compactness_curve = 1

temper_iso_weights = iso_weights ./ 10

betas = iso_weights

ideal_pop = get_node_population(Gct, Gct.root.global_ID) / num_dists
pop_mins_maxs = [(ideal_pop*(1-pop_dev), ideal_pop*(1+pop_dev)) for pop_dev in pop_devs]

constraints = []
measures = []
temper_measures = []
for j=1:num_levels 
	constraint = initialize_constraints() 
	add_constraint!(constraint, PopulationConstraint(pop_mins_maxs[j][1], pop_mins_maxs[j][2]))
	add_constraint!(constraint, ConstrainDiscontinuousTraversals())
	push!(constraints, constraint)

	measure = Measure(0.0)
	push_energy!(measure, compactnessEnergy(compactness_curve), iso_weights[j], "compactness")
	push!(measures, measure)

	temper_measure = Measure(0.0)
	push_energy!(temper_measure, compactnessEnergy(compactness_curve), temper_iso_weights[j], "compactness")
	push!(temper_measures, temper_measure)
end   


rank = 1
processor_to_part_ID = [collect(1:num_levels)]

# ============== OUTPUT SETUP 
writers::Vector{Union{SimpleWriter, Nothing}} = [nothing for _=1:num_levels]
output_file_path_1 = "../output/$(n)x$(n)_districts/"*suffix*".jsonl"
writers[1] = SimpleWriter(measures[1], constraints[1], num_dists, output_file_path_1)

part_ID_to_temp::Dict{Int, Int} = Dict([x => x for x in 1:num_levels])
my_partitions = Dict{Int, MultiLevelPartition}() #part_ID to partition 

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

outer_steps = 1
inner_steps = 10
println("total SNF steps: ", outer_steps*inner_steps*SNF_step_sizes[end], "; total PT flip attempts: ", outer_steps*inner_steps)

total_steps_taken = Dict([x => 0 for x in 1:num_levels])

no_merges = [0 for i=1:num_levels]
no_splits = [0 for i=1:num_levels]
almosts = [0 for i=1:num_levels] 
succs = [0 for i=1:num_levels]
all_ESS_stuff = [[] for _=1:num_levels]
profile_SNF_time = [0.0 for _=1:num_levels] 
profile_swap_time = 0.0

for repeat = 1:outer_steps
	for i=1:inner_steps 
		if mod(i, 500) == 0 
			@show i, now()
			flush(stdout)
		end
	
		for part_ID in 1:num_levels
			t = part_ID_to_temp[part_ID] 
			SNF_step_size = SNF_step_sizes[t]

			for _ =1:SNF_step_size
				isos = []
				for di = 1:num_dists
					push!(isos, get_perimeter(my_partitions[part_ID],di) ^2 / get_area(my_partitions[part_ID].subgraphs[di]))
				end
				sort!(isos)

				lnsp = 0
				for dist_subgraph in my_partitions[part_ID].subgraphs
					s_graph = child_graph(dist_subgraph, dist_subgraph.root.global_ID)[1]
					lnsp += log_nspanning(s_graph)
				end
				push!(isos, lnsp)

				push!(all_ESS_stuff[part_ID], isos)

				tick() 
				run_metropolis_hastings_SNF!(my_partitions[part_ID], measures[t], temper_measures[t],
					(total_steps_taken[part_ID], total_steps_taken[part_ID]+1), constraints[t], rng, writer=writers[t], output_freq=25)
				elapsed = tock() 
				profile_SNF_time[t] += elapsed 
			end

			total_steps_taken[part_ID] += SNF_step_size
			push!(my_partitions[part_ID].extensions["temperature"], t)
		end 
		tick() 
		PT_swap!(TH, my_partitions, part_ID_to_temp, num_levels, pop_mins_maxs, measures, i, rng, betas, 
			(no_merges, no_splits, almosts, succs))
		elapsed = tock() 
		profile_swap_time += elapsed 
		#@show repeat, now()
    end
end

for i=1:num_levels
	if length(all_ESS_stuff[i]) != 0
		open(suffix*"_"*string(i)*"_autocorr_stuff.csv", "w+") do file
			write(file, join([join(z, ",") for z in all_ESS_stuff[i]], "\n"))
		end
	end
end

@show no_merges,no_splits,almosts,succs