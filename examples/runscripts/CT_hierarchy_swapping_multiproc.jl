import Pkg
#codepath = joinpath("/home/postdoc/gtc17/MSMS_copy", "MultiScaleMapSampler") 
codepath = joinpath("../") 
Pkg.activate("mergeSplit"; shared=true)
Pkg.develop(path=codepath)

using RandomNumbers
using Random
using MultiScaleMapSampler
using Test
using JSON 
using CSV 

using TickTock
using Dates

using MPI
using StatProfilerHTML

MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm) + 1
tasks = MPI.Comm_size(comm)


function run_swaps() 
    # ========== GRAPH SETUP 

    CT_json = joinpath("test_graphs", "CT_pct20.json")
    CT_node_data = Set(["COUNTY", "NAME", "POP20", "area", "border_length"])
    CT_base_graph = BaseGraph(CT_json, "POP20", inc_node_data=CT_node_data,
        area_col="area",node_border_col="border_length", edge_perimeter_col="length")
    for k in CT_base_graph.node_attributes 
        k["county_and_pct"] = string(k["COUNTY"], ",",k["NAME"])
    end 

    Gct = MultiLevelGraph(CT_base_graph, ["county_and_pct"]);
    oriented_nbrs = JSON.parsefile("test_graphs/CT_pct20_oriented_nbrs.json")
    oriented_nbrs = Dict{GlobalID, Vector{GlobalID}}([parse(Int, k) => v for (k,v) in oriented_nbrs]);

    filename_prefix = "test_graphs/CTTH_17/tempered_hierarchy_"
    TH = TemperingHierarchy(Gct, filename_prefix, 24);
    num_levels = length(TH.graphs_by_level)
    num_dists = 5

    rngseed = length(ARGS) > 1 ? parse(Int64, ARGS[2]) : 15192
    rng = PCG.PCGStateOneseq(UInt64, rngseed)
    suffix = ARGS[1]
    outer_steps = 400
    inner_steps = 1000

    SNF_step_sizes = [30 for _=1:num_levels]
	
	if rank == 1
		println(filename_prefix)
		println(suffix, " ", outer_steps, " ", inner_steps, " ", rngseed, " ", num_levels)
		println("compactness: -0 -0 -0 -0.5 -0.7")
	end

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

    # ============= MPI setup 
    processor_to_part_ID = [[] for _=1:16]
	for i=1:num_levels 
        push!(processor_to_part_ID[mod(i, 16)+1], i)
    end 

    # ============== OUTPUT SETUP 
    writers::Vector{Union{SimpleWriter, Nothing}} = [nothing for _=1:num_levels]
    output_file_path_1 = "../new_output/"*suffix*string(rank)*"_1.jsonl"
    writers[1] = SimpleWriter(measures[1], constraints[1], num_dists, output_file_path_1)

    my_responsibilities = processor_to_part_ID[rank] 
    part_ID_to_temp::Dict{Int, Int} = Dict([x => x for x in my_responsibilities])
    my_partitions = Dict{Int, MultiLevelPartition}() #part_ID to partition 
    

    for part_ID in my_responsibilities
    	t = part_ID_to_temp[part_ID]
		partition = MultiLevelPartition(TH.graphs_by_level[t], constraints[t], num_dists; rng=rng, max_attempts=100)
		populate_extras!(TH.graphs_by_level[t])
        partition.extensions["iso"] = []; partition.extensions["aps"] = []; partition.extensions["pop"] = []
        partition.extensions["temperature"] = []
        partition.extensions["partition_ID"] = part_ID
        my_partitions[part_ID] = partition

        run_metropolis_hastings_SNF!(my_partitions[part_ID], measures[t], temper_measures[t], 
            1000, constraints[t], rng, writer=writers[t], output_freq=100)
    end
    MPI.Barrier(comm)

	total_steps_taken = Dict([x => 0 for x in my_responsibilities])

    if rank == 1 
		tick() 
	end

	no_merges = [0 for i=1:num_levels]
	no_splits = [0 for i=1:num_levels]
	almosts = [0 for i=1:num_levels] 
	succs = [0 for i=1:num_levels]
	NSNs = [[] for i=1:num_levels]
	ps = [[] for i=1:num_levels]
	pps = [([],[],[],[]) for i=1:num_levels]	

    for repeat = 1:outer_steps
        for i=1:inner_steps 
            if mod(i, 200) == 0 
                for part_ID in my_responsibilities 
                	t = part_ID_to_temp[part_ID]
                	if t == 1
                    	write_TH_partition_to_file(TH, my_partitions[part_ID], t, 
                        	"viz/tempering_hierarchy2/"*suffix*string(t)*"partition_"*string(total_steps_taken[part_ID])*".csv")
               		end
                end
                if rank == 1
                    @show i, now()
                    flush(stdout)
                end
            end

            for part_ID in my_responsibilities
                t = part_ID_to_temp[part_ID] 
                SNF_step_size = SNF_step_sizes[t]
                run_metropolis_hastings_SNF!(my_partitions[part_ID], measures[t], temper_measures[t], 
                    (total_steps_taken[part_ID], total_steps_taken[part_ID]+SNF_step_size), constraints[t], rng, writer=writers[t], output_freq=25)
                total_steps_taken[part_ID] += SNF_step_size
                push!(my_partitions[part_ID].extensions["temperature"], t)
            end 

            # (1). Set finers and coarsers 
            finers   = mod(i, 2) == 0 ? (1:2:(num_levels-1)) : (2:2:(num_levels-1)) 
            coarsers = mod(i, 2) == 0 ? (2:2:num_levels) : (3:2:num_levels)

            # (2). compute split node counts of each partition in a fnr level 
            split_nodes_by_level::Vector{Int} = [0 for _ = 1:(num_levels-1)]
            for part_ID in my_responsibilities 
                t = part_ID_to_temp[part_ID]
                if t in finers 
                    split_nodes_by_level[t] = count_split_nodes(my_partitions[part_ID], TH.children_by_level[t+1])
                end
            end
            # and allreduce them so everyone knows all of them. 
            MPI.Allreduce!(MPI.IN_PLACE, split_nodes_by_level, MPI.SUM, comm)

            potential_swaps = Dict{Int, MultiLevelPartition}()
            # (3) compute the swap probabilities and stuff! 
            log_p1_factors = [0.0 for _ = 1:num_levels]
            log_p2_factors = [0.0 for _ = 1:num_levels]
            for part_ID in my_responsibilities 
                t = part_ID_to_temp[part_ID]
                if t in finers 
                    fnr, csr = t, t+1
                    f2c_partition = copy_partition(my_partitions[part_ID])
                    log_p_up = swap_up(f2c_partition, TH, fnr, pop_mins_maxs, split_nodes_by_level[fnr], rng, betas[fnr])

                    if log_p_up === nothing
                        log_p1_factors[fnr] = -Inf 
                        log_p2_factors[fnr] = -Inf 
                        no_merges[fnr] += 1 
                    else 
                        cpartition_new = project_up_pctpair(f2c_partition, TH.graphs_by_level[csr], TH.parents_by_level[fnr]) 
                        cur_energy = fully_compute_energy(my_partitions[part_ID], measures[fnr]) 
                        new_energy = fully_compute_energy(cpartition_new, measures[csr])
                        log_p1_factors[fnr] = log_p_up
                        log_p2_factors[fnr] = (-new_energy + cur_energy)
                        potential_swaps[part_ID] = cpartition_new
                    end
                elseif t in coarsers
                    fnr, csr = t-1, t
                    c2f_partition = project_down_pctpair(my_partitions[part_ID], 
                        TH.graphs_by_level[fnr], TH.children_by_level[csr])

                    log_p_down = swap_down(c2f_partition, TH, fnr, pop_mins_maxs, split_nodes_by_level[fnr], rng, betas[fnr]) 

                    if log_p_down === nothing 
                        log_p1_factors[csr] = -Inf
                        log_p2_factors[csr] = -Inf
                        no_splits[csr] += 1
                    else 
                        cur_energy = fully_compute_energy(my_partitions[part_ID], measures[csr]) 
                        new_energy = fully_compute_energy(c2f_partition, measures[fnr])
                        log_p1_factors[csr] = log_p_down
                        log_p2_factors[csr] = (-new_energy + cur_energy)
                        potential_swaps[part_ID] = c2f_partition
                    end
                else 
                    # one of the ends, and we're not swapping it this time 
                end 
            end

            # (4) combine the swap probabilities, and decide which temps need to swap 
            MPI.Allreduce!(MPI.IN_PLACE, log_p1_factors, MPI.SUM, comm)
            MPI.Allreduce!(MPI.IN_PLACE, log_p2_factors, MPI.SUM, comm)
            log_p_factors = log_p1_factors .+ log_p2_factors 

            if false # rank == 1 
                for t in finers
                    p = exp(log_p_factors[t] + log_p_factors[t+1])
                    push!(ps[t], p)
                    push!(pps[t][1], log_p1_factors[t])
                    push!(pps[t][2], log_p2_factors[t])
                    push!(pps[t][3], log_p1_factors[t+1])
                    push!(pps[t][4], log_p2_factors[t+1])

                end
            end

            # the finer indices will be true if the swap happens 
            swaps = [false for _=1:num_levels]
            for part_ID in my_responsibilities 
                t = part_ID_to_temp[part_ID]
                if t in finers
                    p = exp(log_p_factors[t] + log_p_factors[t+1])
                    if !isapprox(p, 0.0) 
                        almosts[t] += 1
                    end
                    if rand(rng) < p
                        succs[t] += 1
                        swaps[t] = true 
                    end
                end
            end
            MPI.Allreduce!(MPI.IN_PLACE, swaps, MPI.LOR, comm)

            for part_ID in my_responsibilities 
                t = part_ID_to_temp[part_ID]
                if t in finers && swaps[t]
                    part_ID_to_temp[part_ID] = t+1
                    my_partitions[part_ID] = potential_swaps[part_ID]

                    run_metropolis_hastings_SNF!(my_partitions[part_ID], measures[t+1], temper_measures[t+1], 
                        (total_steps_taken[part_ID], total_steps_taken[part_ID]+50), constraints[t+1], rng, writer=nothing, output_freq=25)
                    total_steps_taken[part_ID] += 50
                elseif t in coarsers && swaps[t-1]
                    part_ID_to_temp[part_ID] = t-1
                    my_partitions[part_ID] = potential_swaps[part_ID]

                    run_metropolis_hastings_SNF!(my_partitions[part_ID], measures[t-1], temper_measures[t-1], 
                        (total_steps_taken[part_ID], total_steps_taken[part_ID]+50), constraints[t-1], rng, writer=nothing, output_freq=25)
                    total_steps_taken[part_ID] += 50 
                end
            end
        end
        if rank == 1
            @show repeat, now()
            #@show d
        end	
        for (part_ID, part) in my_partitions
            temps_s = join(part.extensions["temperature"], ",") * ","
            part.extensions["temperature"] = []
            open(joinpath("..", "new_output", string(part_ID)*suffix*"temperatures_MP.csv"), "a") do file 
                write(file, temps_s)
            end
        end

    end
    if rank == 1 
        #@show SNBLs
        tock() 
    end


    for (part_ID, part) in my_partitions
        t = part_ID_to_temp[part_ID]
        #@show part.extensions
        open(joinpath("..", "new_output", string(part_ID)*suffix*"temperatures_MP.csv"), "a") do file 
            write(file, "\n")
        end

        write_partition_to_file_flat(my_partitions[part_ID],
            "../new_output/"*suffix*string(part_ID)*"_final_partition.csv")
        
		open(joinpath("..", "new_output", string(part_ID)*suffix*"extensions.json"), "a") do file 
            JSON.print(file, part.extensions)
        end
    end

	MPI.Allreduce!(MPI.IN_PLACE, no_splits, MPI.SUM, comm) 
	MPI.Allreduce!(MPI.IN_PLACE, no_merges, MPI.SUM, comm) 
	MPI.Allreduce!(MPI.IN_PLACE, almosts, MPI.SUM, comm) 
	MPI.Allreduce!(MPI.IN_PLACE, succs, MPI.SUM, comm) 
	
	open(joinpath("../new_output/"*suffix*"_swapinfo.csv"), "w+") do file 
		write(file, join(no_splits, ",") * "\n")
		write(file, join(no_merges, ",") * "\n")
		write(file, join(almosts, ",") * "\n")
		write(file, join(succs, ",") * "\n")
	end

    println(rank, " is done.")
end 

run_swaps() 
