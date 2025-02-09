""""""
function run_metropolis_hastings!(
    partition::MultiLevelPartition,
    proposal::Union{Function,Vector{Tuple{T, Function}}},
    measure::Measure,
    steps::Union{Int,Tuple{Int,Int}},
    rng::AbstractRNG;
    writer::Union{Nothing, SimpleWriter}=nothing, 
    output_freq=25
) where T <: Real
    #unclear if this is needed 
    #precompute_node_tree_counts!(partition)

    #if the proposal is a vector, check that the probabilities sum to 1 
    if typeof(proposal) <: Vector 
        @assert sum(proposal[i][1] for i = 1:length(proposal)) == 1 
    end

    initial_step, final_step = typeof(steps)<:Tuple ? steps : (1, steps)
    if initial_step == 0 || initial_step == 1
        output(partition, measure, initial_step, 0, writer)
    end

    step = initial_step 
    while step <= final_step 
        proposal!, proposal_index = get_random_proposal(proposal, rng)

        p, update = proposal!(partition, measure, rng, step)
        # update is a tuple: (
            # (d1::Int, d2::Int), 
            # new_dist_masks (mask1, mask2) 
            # new_dist_pops (pop1, pop2) 
            # edge (UndirectedEdge(GlobalID, GlobalID))
        #)

        # for clarity: The `p` that returned by the proposal should be generally considered 
        # to be the ratio of the forward and backward proposal probabilities. 
        # That is, when metropolizing, we should in general take this `p` and multiply it by the 
        # pi(proposed state) / pi(current state) ratio. 
        # however, when gamma is less than 1, our overall measure includes some tree count 
        # and linking edge count nonsense. In these cases, this p already includes that 
        # part of the measure. 

        push!(partition.extensions[TRANSITION_RATIOS], p)

        if update === nothing || p == 0 
            #println("update is nothing")
            if mod(step, output_freq) == 0 && step != initial_step
                output(partition, measure, step, 0, writer)
            end
            partition.extensions[REJECTION_COUNT]["rejection_p0"]+=1
            step += 1
            continue
        end
  
        de = get_delta_energy(partition, measure, update)
        p *= de
        if rand(rng) < p
            #println("accepting update, p=", p) 
            partition.extensions[ACCEPTANCE_COUNT] += 1
            if typeof(update) <: Tuple # forest recom update
                update_partition!(partition, update)
            elseif typeof(update) <: Function # single node flip update 

            end
        else
            #println("rejecting update, p=", p)
            partition.extensions[REJECTION_COUNT]["rejection_metropolis"]+=1
        end
        
        if mod(step, output_freq) == 0 && step != initial_step
            output(partition, measure, step, 0, writer)
        end
        step += 1 
    end
end

function run_metropolis_hastings_SNF!(
    partition::MultiLevelPartition,
    measure::Measure,
    temper_measure::Measure,
    steps::Union{Int,Tuple{Int,Int}},
    constraints, 
    rng::AbstractRNG;
    writer::Union{Nothing, SimpleWriter}=nothing, 
    output_freq=25, 
    output_freq_offset=0,
    print=false
)

    #measure = Measure(0.0)

    initial_step, final_step = typeof(steps)<:Tuple ? steps : (1, steps)
    if initial_step == 0 || initial_step == 1
        output(partition, measure, initial_step, 0, writer)
    end

    step = initial_step 
    while step <= final_step 
        p_fwd, reject_update, accept_update, _, pct, dist_to = 
            single_node_flip_basic(partition, temper_measure, constraints, rng)        
        dist_from = look_up_district(partition, pct)

        old_aps = Dict{Int, Set{GlobalID}}([
            dist_from => partition.articulation_points[dist_from], 
            dist_to => partition.articulation_points[dist_to]])

        accept_SNF!(partition, dist_to, pct)
        new_energy = fully_compute_energy(partition, measure)
        p_bck = single_node_flip_basic(partition, temper_measure, constraints, rng, just_prob_of_flip=(dist_from, pct))

        new_aps = Dict{Int, Set{GlobalID}}([
            dist_from => partition.articulation_points[dist_from], 
            dist_to => partition.articulation_points[dist_to]])

        accept_SNF!(partition, dist_from, pct; known_artpts=old_aps)

        cur_energy = fully_compute_energy(partition, measure)
        #@show new_energy, cur_energy, pct, dist_to, dist_from, p_bck, p_fwd 
        p = (p_bck / p_fwd) * exp(-new_energy + cur_energy)
        #a6 += @allocated push!(partition.extensions[TRANSITION_RATIOS], p_bck / p_fwd)
        #a6 += @allocated push!(partition.extensions[TRANSITION_RATIOS2], exp(-new_energy + cur_energy))

        if print
            @show 1/p_fwd, 1/p_bck, cur_energy, new_energy, p
        end

        if rand(rng) < p
            if print println("accepting update, p=", p) end
            partition.extensions[ACCEPTANCE_COUNT] += 1
            accept_update(new_aps) 
        else
            if print println("rejecting update, p=", p) end
            partition.extensions[REJECTION_COUNT]["rejection_metropolis"]+=1
            reject_update(nothing)
        end
        
        if mod(step+output_freq_offset, output_freq) == 0 && step != initial_step
            output(partition, measure, step, 0, writer)
        end
        step += 1 

        #@assert all([length(mask[partition.graph.root.global_ID]) > 1 for mask in partition.district_to_nodes])
    end
end

function get_random_proposal(
    proposal::Union{Function, Vector{Tuple{T, Function}}},
    rng::AbstractRNG
) where T <: Real
    if typeof(proposal) <: Function 
        return proposal, 1
    end
    proposal_weights = [i[1] for i in proposal]
    index = findfirst(cumsum(proposal_weights) .> rand(rng))
    return proposal[index][2], index
end

function update_partition!(
    partition::MultiLevelPartition,
    update::Tuple
)::Nothing
    changed_districts, new_dist_masks, new_dist_pops, edge = update
    old_n2d = length(partition.node_to_district)

    clear_node_to_district!(partition, changed_districts)
    clear_adjacency!(partition, changed_districts)
    clear_shared_nodes!(partition, changed_districts)
    clear_linking_edge_attrs!(partition, changed_districts)

    for (ii,cd) in enumerate(changed_districts)
        partition.district_to_nodes[cd] = new_dist_masks[ii]
        partition.subgraphs[cd] = take_subgraph(partition.graph, new_dist_masks[ii])
        partition.dist_populations[cd] = new_dist_pops[ii]
    end
    set_node_to_district!(partition, changed_districts)
    set_adjacency_and_shared_nodes!(partition)
    set_linking_edge_attrs!(partition)

    for (ii,cd) in enumerate(changed_districts)
        partition.dist_perimeters[cd] = get_perimeter(partition, cd)
        partition.dist_areas[cd] = get_area(partition.subgraphs[cd])
    end
end



function output(
    partition::MultiLevelPartition,
    measure::Measure,
    step::Int,
    count::Int,
    writer::Union{SimpleWriter, Nothing}
)
    if writer === nothing
        return
    end

    if typeof(writer) <: SimpleWriter 
        d_initial = Dict{String, Any}("step" => step, "weight" => 1, "name" => "step"*string(step), "data" => Dict())
        write_districts(writer, partition, d_initial)

    else
        @assert false 
        
        for (desc, f) in writer.map_output_data
            writer.map_param[desc] = f(partition)
        end
        if !writer.output_districting
            d = Dict{Tuple{Vararg{String}}, Int}()
            map = Map{MapParam}("step"*string(step-count), d, 1, writer.map_param)
        else
            # partition.node_to_district is a dict from GlobalID to district #. 
            # for Map, we want it to be a dict from variadic-length tuples of strings to district #. 
            # for now, it's a vector (cuz who uses variadic length tuples anyway???) 
            n2d_for_output = Dict{Tuple{Vararg{String}},Int64}() 

            for dist in 1:partition.num_dists 
                function f(G, k)
                    node_name = node_name_exclude_quotient_nodes(G, k.global_ID)
                    if k.node_name !== nothing 
                        if length(node_name) == 1 && occursin(",", node_name[1])
                            node_name = split(node_name[1], ",")
                        end 
                        name_tuple = Tuple{Vararg{String}}(node_name)
                        n2d_for_output[name_tuple] = dist
                    end
                end
            
                for (k,v) in partition.district_to_nodes[dist] 
                    if v === nothing 
                        if partition.graph.ID_to_node[k].node_name !== nothing 
                            node_name = node_name_exclude_quotient_nodes(partition.graph, k)
                            if length(node_name) == 1 && occursin(",", node_name[1])
                                node_name = split(node_name[1], ",")
                            end 
                            name_tuple = Tuple{Vararg{String}}(node_name)
                            n2d_for_output[name_tuple] = dist    
                        elseif length(node_name_exclude_quotient_nodes(partition.graph, k)) <= 1
                            preord_foreach(partition.subgraphs[dist], k, f)
                        end
                    end
                end
            end

            map = Map{MapParam}("step"*string(step-count), 
                n2d_for_output, 1, writer.map_param)
        end

        try
            addMap(writer.atlas.io, map)
        catch e
            println("Could not add map to atlas")
            @assert false
        end
        if haskey(partition.extensions, rejection_counter::EXTENSIONS)
            partition.extensions[rejection_counter::EXTENSIONS]["acceptance_wait"] = 0
        end
    end 
end 
