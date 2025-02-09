function pop_balance(dist_pops, min, max) 
    total = 0 
    for dpop in dist_pops 
        if min <= dpop && dpop <= max continue 
        elseif min > dpop 
            total += min - dpop 
        else 
            total += dpop - max 
        end
    end
    return total 
end

function all_dist_art_pts(partition) 
    art_pts = Set{GlobalID}()
    for (_,x) in partition.articulation_points
        union!(art_pts, x)
    end
    return art_pts 
end 

function project_up_pctpair(partition, Gncpair, pair_parents)
    root = partition.graph.root.global_ID
    new_d2n = []
    for mask in partition.district_to_nodes 
        m2 = Dict{GlobalID, Union{Nothing, Set{GlobalID}}}() 
        root_children = Set{GlobalID}()
        for node in mask[root]
            parent = pair_parents[node]
            m2[parent] = nothing 
            push!(root_children, parent)
        end
        m2[Gncpair.root.global_ID] = root_children
        push!(new_d2n, m2)
    end
    return partition_from_d2n(Gncpair, new_d2n, partition.num_dists; extensions=partition.extensions, care_about_linking_edges=false)
end

function project_down_pctpair(partition::MultiLevelPartition, Gncpct, pair_children) 
    root = partition.graph.root.global_ID
    new_d2n = []
    for mask in partition.district_to_nodes 
        m2 = Dict{GlobalID, Union{Nothing, Set{GlobalID}}}() 
        root_children = Set{GlobalID}()
        for node in mask[root]
            for child in pair_children[node]
                m2[child] = nothing 
                push!(root_children, child)
            end
        end
        m2[Gncpct.root.global_ID] = root_children
        push!(new_d2n, m2)
    end
    return partition_from_d2n(Gncpct, new_d2n, partition.num_dists; extensions=partition.extensions, care_about_linking_edges=false)
end

function count_split_nodes(partition, pair_children) 
    ct = 0
    for (_, pair) in pair_children 
        c = Set([look_up_district(partition, x) for x in pair])
        if length(c) > 1 
            ct += 1 
        end 
    end
    return ct
end

function get_valid_splits_pair_tempered(
    Gncpct,
    fine_partition, 
    fine_min, fine_max,
    coarse_min, coarse_max, 
    balance_limit,
    pair_parents, pair_children, 
    iso_limit, rng, beta; look_for = nothing
)
    valid_moves = Dict{Tuple{GlobalID, Int, Int}, Real}()
    total_weight = 0
    sibling_node(z) = pair_children[pair_parents[z]][1] == z ? pair_children[pair_parents[z]][2] : pair_children[pair_parents[z]][1]
    is_singleton(z) = length(pair_children[pair_parents[z]]) == 1

    art_pts = all_dist_art_pts(fine_partition)

    og_iso = sum((fine_partition.dist_perimeters .^ 2) ./ fine_partition.dist_areas)

    for node_ID in finest_level_nodes(Gncpct; is_flat=true)
        d1 = look_up_district(fine_partition, node_ID)
        if is_singleton(node_ID) 
            continue # only child of its parent 
        end
        if look_up_district(fine_partition, sibling_node(node_ID)) != d1
            continue # this node's parent coarse node is already split 
        end 
        if node_ID in art_pts
            continue #this node is an articulation point and can't be moved 
        end 

        for d2 in dists_adj_to_node(fine_partition, node_ID)
            npop = get_node_population(Gncpct, node_ID)
            fine_partition.dist_populations[d1] -= npop 
            fine_partition.dist_populations[d2] += npop 
            new_fine_balance = pop_balance(fine_partition.dist_populations, fine_min, fine_max) 
            new_coarse_balance = pop_balance(fine_partition.dist_populations, coarse_min, coarse_max) 

            if new_fine_balance <= balance_limit && new_coarse_balance == 0
                flip_iso = iso_if_you_flipped(fine_partition, node_ID, d1, d2)
                if iso_limit === nothing || flip_iso <= iso_limit    
                    weight = exp(-beta*(flip_iso - og_iso))
                    total_weight += weight 
                    valid_moves[node_ID, d1, d2] = weight
                end
            end
            fine_partition.dist_populations[d1] += npop 
            fine_partition.dist_populations[d2] -= npop 
        end
    end

    if length(valid_moves) == 0 
        return nothing 
    end
    if look_for !== nothing 
        return valid_moves[look_for] / total_weight 
    end

    choice = rand(rng) * total_weight
    so_far = 0 
    for (k, weight) in valid_moves 
        node_ID, d1, d2 = k
        so_far += weight 
        if so_far > choice
            return node_ID, d1, d2, weight/total_weight 
        end
    end
    @assert false 
end

function iso_if_you_flipped(partition, pct, dist_from, dist_to)
    merged_nodes = Dict{GlobalID, Set{GlobalID}}()
    pct_perim_dist_from = node_perim_with_district(partition, merged_nodes, pct, dist_from)
    pct_perim_dist_to = node_perim_with_district(partition, merged_nodes, pct, dist_to)
    pct_perim_total = node_perim(partition, pct; merged_nodes=merged_nodes)
    pct_area = get_node_area(partition.graph, pct; merged_nodes=merged_nodes) 

    areas, perims = partition.dist_areas, partition.dist_perimeters 
    perims[dist_from] -= (pct_perim_total - 2*pct_perim_dist_from)
    perims[dist_to] += (pct_perim_total - 2*pct_perim_dist_to)
    areas[dist_from] -= pct_area 
    areas[dist_to] += pct_area 
    
    iso = 0

    if !haskey(partition.extensions, "norm")
        iso = sum(perims .^ 2 ./ areas)
    else 
        norm = partition.extensions["norm"]
        iso = sum((perims .^ 2 ./ areas) .^ norm ) .^ (1/norm)
    end
        
    perims[dist_from] += (pct_perim_total - 2*pct_perim_dist_from)
    perims[dist_to] -= (pct_perim_total - 2*pct_perim_dist_to)
    areas[dist_from] += pct_area 
    areas[dist_to] -= pct_area 

    return iso 
end    

function get_valid_merges_pair(
    Gncpct, Gncpair,
    partition, 
    fine_min, fine_max, 
    coarse_min, coarse_max, 
    balance_limit, 
    pair_parents, pair_children
)
    art_pts = all_dist_art_pts(partition) 

    function is_one_dist(z) 
        if length(pair_children[z]) == 1 
            return true 
        else 
            return look_up_district(partition, first(pair_children[z])) == look_up_district(partition, last(pair_children[z]))
        end
    end

    #child_dists(z) = Set([look_up_district(partition, x) for x in pair_children[z]])
    valid_moves = [] 
    for node_ID in finest_level_nodes(Gncpair; is_flat=true)
        if is_one_dist(node_ID) #length(child_dists(node_ID)) == 1 
            continue 
        end 

        n1, n2 = pair_children[node_ID]
        d1, d2 = look_up_district(partition, n1), look_up_district(partition, n2)
    
        this_node_valid_moves = []

        if !(n1 in art_pts) 
            # flip n1 to d2, that is, make node_ID wholly d2 

            npop = get_node_population(Gncpct, n1)
            partition.dist_populations[d1] -= npop
            partition.dist_populations[d2] += npop
            new_fine_balance = pop_balance(partition.dist_populations, fine_min, fine_max) 
            new_coarse_balance = pop_balance(partition.dist_populations, coarse_min, coarse_max) 
            partition.dist_populations[d1] += npop
            partition.dist_populations[d2] -= npop
            if new_fine_balance <= balance_limit && new_coarse_balance == 0
                push!(this_node_valid_moves, (n1, d1, d2, iso_if_you_flipped(partition, n1, d1, d2) ))
            end
        end    
        if !(n2 in art_pts) 
            # flip n2 to d1, that is, make node_ID wholly d1
            npop = get_node_population(Gncpct, n2) 
            partition.dist_populations[d1] += npop
            partition.dist_populations[d2] -= npop
            new_fine_balance = pop_balance(partition.dist_populations, fine_min, fine_max) 
            new_coarse_balance = pop_balance(partition.dist_populations, coarse_min, coarse_max) 
            partition.dist_populations[d1] -= npop
            partition.dist_populations[d2] += npop

            if new_fine_balance <= balance_limit && new_coarse_balance == 0
                push!(this_node_valid_moves, (n2, d2, d1, iso_if_you_flipped(partition, n2, d2, d1)))
            end
        end

        # when a given node can be moved either way, the only one we are allowed to use is the 
        # one with a better iso score. Is this okay? 
        # this means that when we split a node we're only allowed to split if merging back is better 
        # than merging "forward". Which is Not Necessarily a case where we will screw our 
        # forward iso, but it might be. 

        # This is maybe a bit concerning? I don't know...
        if length(this_node_valid_moves) == 1 
            push!(valid_moves, (this_node_valid_moves[1][1:3], 1))
            #push!(valid_moves, (this_node_valid_moves[1][1:3], 4))
            #push!(valid_moves, (this_node_valid_moves[1][1:3], 4))
            #push!(valid_moves, (this_node_valid_moves[1][1:3], 4))
        elseif length(this_node_valid_moves) == 2 
            if this_node_valid_moves[1][4] < this_node_valid_moves[2][4]
                #push!(valid_moves, (this_node_valid_moves[1][1:3], 3))
                #push!(valid_moves, (this_node_valid_moves[1][1:3], 3))
                push!(valid_moves, (this_node_valid_moves[1][1:3], 1))
                push!(valid_moves, (this_node_valid_moves[2][1:3], 1))
            else
                push!(valid_moves, (this_node_valid_moves[1][1:3], 1))
                push!(valid_moves, (this_node_valid_moves[2][1:3], 1))
                #push!(valid_moves, (this_node_valid_moves[2][1:3], 3))
                #push!(valid_moves, (this_node_valid_moves[2][1:3], 3))
            end
        end
    end
    return valid_moves 
end

#returns log_p (total backward prob / total forward prob)
function swap_up(f2c_partition::MultiLevelPartition, TH::TemperingHierarchy, fnr::Int, 
    pop_mins_maxs, num_steps_to_take::Int, rng, beta
)
    csr = fnr + 1
    fine_min, fine_max = pop_mins_maxs[fnr]
    coarse_min, coarse_max = pop_mins_maxs[csr]
    pair_children = TH.children_by_level[csr]
    pair_parents = TH.parents_by_level[fnr]

    G_fnr = TH.graphs_by_level[fnr]
    G_csr = TH.graphs_by_level[csr]

    max_pop_unbalance = f2c_partition.num_dists * ((coarse_max - fine_max) + (fine_min - coarse_min))/2
    min_pop_unbalance = 0.0
    balance_limit_step_size = (max_pop_unbalance - min_pop_unbalance) / num_steps_to_take
    f2c_balance_limit = min_pop_unbalance

    #beta = 0.3
    log_p = 0

    for _ = 1:num_steps_to_take
        @assert pop_balance(f2c_partition.dist_populations, fine_min, fine_max) <= f2c_balance_limit
        @assert pop_balance(f2c_partition.dist_populations, coarse_min, coarse_max) == 0
      
        f2c_balance_limit += balance_limit_step_size
        #@show f2c_balance_limit
        merges = get_valid_merges_pair(G_fnr, G_csr, f2c_partition, fine_min, fine_max, coarse_min, coarse_max, 
            f2c_balance_limit, pair_parents, pair_children)

        if length(merges) == 0 
            return nothing # no valid merges 
        end

        # pick a uniform merge 
        merge_to_make, p_choose_this_merge = uniformly_choose(merges, rng) 
        ((merge_pct, merge_dist_from, merge_dist_to), multiplier) = merge_to_make 
        p_choose_this_merge *= multiplier
        
        # make the moves! 
        accept_SNF!(f2c_partition, merge_dist_to, merge_pct)

        # backward probability 
        p_choose_this_bwd_split = get_valid_splits_pair_tempered(G_fnr, f2c_partition, fine_min, fine_max, 
            coarse_min, coarse_max, f2c_balance_limit - balance_limit_step_size, 
            pair_parents, pair_children, nothing, rng, beta; look_for=(merge_pct, merge_dist_to, merge_dist_from))

        p_incremental = p_choose_this_bwd_split / p_choose_this_merge
        log_p += log(p_incremental)
    end
    @assert pop_balance(f2c_partition.dist_populations, coarse_min, coarse_max) == 0

    return log_p  
end

function swap_down(c2f_partition::MultiLevelPartition, TH::TemperingHierarchy, fnr::Int, 
    pop_mins_maxs, num_steps_to_take, rng, beta
)
    csr = fnr + 1
    fine_min, fine_max = pop_mins_maxs[fnr]
    coarse_min, coarse_max = pop_mins_maxs[csr]
    pair_children = TH.children_by_level[csr]
    pair_parents = TH.parents_by_level[fnr]

    G_fnr = TH.graphs_by_level[fnr]
    G_csr = TH.graphs_by_level[csr]

    max_pop_unbalance = c2f_partition.num_dists * ((coarse_max - fine_max) + (fine_min - coarse_min))/2
    min_pop_unbalance = 0.0
    balance_limit_step_size = (max_pop_unbalance - min_pop_unbalance) / num_steps_to_take
    c2f_balance_limit = max_pop_unbalance

    #beta = 0.3
    log_p = 0

    if num_steps_to_take == 0 && pop_balance(c2f_partition.dist_populations, fine_min, fine_max) != 0
        return nothing # there are no splits to fix but the pop balance doesn't work 
    end

    for _ = 1:num_steps_to_take
        @assert pop_balance(c2f_partition.dist_populations, fine_min, fine_max) <= c2f_balance_limit
        @assert pop_balance(c2f_partition.dist_populations, coarse_min, coarse_max) == 0
      
        c2f_balance_limit -= balance_limit_step_size 
        #@show c2f_balance_limit
        if (c2f_balance_limit < 0) c2f_balance_limit = 0 end 

        splits = get_valid_splits_pair_tempered(G_fnr, c2f_partition, fine_min, fine_max, coarse_min, coarse_max, 
            c2f_balance_limit, pair_parents, pair_children, nothing, rng, beta)
        if splits === nothing 
            return nothing # no valid splits 
        end
        
        # make the move! 
        split_pct, split_dist_from, split_dist_to, p_choose_this_split = splits
        accept_SNF!(c2f_partition, split_dist_to, split_pct)

        bwd_merges = get_valid_merges_pair(G_fnr, G_csr, c2f_partition, fine_min, fine_max, 
            coarse_min, coarse_max, c2f_balance_limit + balance_limit_step_size, 
            pair_parents, pair_children)
        
        bwd_merge_idx = findfirst(x -> x[1] == (split_pct, split_dist_to, split_dist_from), bwd_merges)
        @assert bwd_merge_idx !== nothing 
        p_choose_this_bwd_merge = bwd_merges[bwd_merge_idx][2]/length(bwd_merges)

        p_incremental = p_choose_this_bwd_merge / p_choose_this_split
        log_p += log(p_incremental)
    end
    @assert pop_balance(c2f_partition.dist_populations, fine_min, fine_max) == 0

    return log_p
end

function swap_pair_pct_tempered(f2c_partition, c2f_partition, Gncpair, Gncpct, rng, 
    fine_min, fine_max, coarse_min, coarse_max, fine_labels, pair_parents, pair_children) 
    no_merged_nodes = Dict{GlobalID, Set{GlobalID}}()

    num_steps_to_take = count_split_nodes(f2c_partition, pair_children)
    max_pop_unbalance = f2c_partition.num_dists * ((coarse_max - fine_max) + (fine_min - coarse_min))/2
    min_pop_unbalance = 0.0
    balance_limit_step_size = (max_pop_unbalance - min_pop_unbalance) / num_steps_to_take

    c2f_balance_limit = max_pop_unbalance
    f2c_balance_limit = min_pop_unbalance

    c_iso = nothing; f_iso = nothing

    beta = 0.3

    p = 1

    if num_steps_to_take == 0 && pop_balance(c2f_partition.dist_populations, fine_min, fine_max) != 0
        return nothing, 3 # there are no splits to fix but the pop balance doesn't work 
    end

    for _ = 1:num_steps_to_take
        @assert pop_balance(f2c_partition.dist_populations, fine_min, fine_max) <= f2c_balance_limit
        @assert pop_balance(c2f_partition.dist_populations, fine_min, fine_max) <= c2f_balance_limit
        @assert pop_balance(f2c_partition.dist_populations, coarse_min, coarse_max) == 0
        @assert pop_balance(c2f_partition.dist_populations, coarse_min, coarse_max) == 0
      
        c2f_balance_limit -= balance_limit_step_size 
        if (c2f_balance_limit < 0) c2f_balance_limit = 0 end 
        f2c_balance_limit += balance_limit_step_size

        splits = get_valid_splits_pair_tempered(Gncpct, c2f_partition, fine_min, fine_max, coarse_min, coarse_max, 
            c2f_balance_limit, pair_parents, pair_children, c_iso, rng, beta)
        merges = get_valid_merges_pair(Gncpct, Gncpair, f2c_partition, fine_min, fine_max, coarse_min, coarse_max, 
            f2c_balance_limit, pair_parents, pair_children)

        if splits === nothing 
            return nothing, 1 # no valid splits 
        elseif length(merges) == 0 
            return nothing, 2, f2c_balance_limit # no valid merges 
        end

        # pick a uniform merge 
        merge_to_make, p_choose_this_merge = uniformly_choose(merges, rng) 
        ((merge_pct, merge_dist_from, merge_dist_to), multiplier) = merge_to_make 
        p_choose_this_merge *= multiplier

        split_pct, split_dist_from, split_dist_to, p_choose_this_split = splits
        
        # make the moves! 
        accept_SNF!(f2c_partition, merge_dist_to, merge_pct)
        accept_SNF!(c2f_partition, split_dist_to, split_pct)

        p_choose_this_bwd_split = get_valid_splits_pair_tempered(Gncpct, f2c_partition, fine_min, fine_max, 
            coarse_min, coarse_max, f2c_balance_limit - balance_limit_step_size, 
            pair_parents, pair_children, f_iso, rng, beta; look_for=(merge_pct, merge_dist_to, merge_dist_from))
        bwd_merges = get_valid_merges_pair(Gncpct, Gncpair, c2f_partition, fine_min, fine_max, 
            coarse_min, coarse_max, c2f_balance_limit + balance_limit_step_size, 
            pair_parents, pair_children)
        
        bwd_merge_idx = findfirst(x -> x[1] == (split_pct, split_dist_to, split_dist_from), bwd_merges)
        @assert bwd_merge_idx !== nothing 
        p_choose_this_bwd_merge = bwd_merges[bwd_merge_idx][2]/length(bwd_merges)

        p_incremental = p_choose_this_bwd_merge * p_choose_this_bwd_split / (p_choose_this_split * p_choose_this_merge)
        p *= p_incremental 
    end
    @assert pop_balance(f2c_partition.dist_populations, coarse_min, coarse_max) == 0
    @assert pop_balance(c2f_partition.dist_populations, fine_min, fine_max) == 0

    return p, num_steps_to_take 
end

function copy_partition(partition; extensions=partition.extensions)
    subgraphs = [take_subgraph(partition.graph, x) for x in partition.district_to_nodes]
    return MultiLevelPartition(
        partition.num_dists, 
        deepcopy(partition.district_to_nodes), 
        deepcopy(partition.node_to_district),
        partition.graph, 
        subgraphs, 
        deepcopy(partition.dist_populations), 
        deepcopy(partition.dist_areas), 
        deepcopy(partition.dist_perimeters), 
        deepcopy(partition.adjacency), 
        deepcopy(partition.shared_node), 
        deepcopy(partition.adjacent_fine_neighbors), 
        deepcopy(partition.linking_edge_attrs), 
        partition.parent, 
        deepcopy(extensions),
        deepcopy(partition.proposed_extensions),
        deepcopy(partition.favorite_pcts), 
        deepcopy(partition.articulation_points)
    )
end

function print_partition(partition, TH) 
    if length(partition.node_to_district) == 17 
        for i=0:3 
            for j=0:3 
                print(partition.node_to_district[i*4+j+1])
            end 
            println() 
        end 
    else 
        @assert length(partition.node_to_district) == 13 
        for i=0:3 
            for j=0:3 
            #if haskey(partition.node_to_district, i*4+j+1)
            #    print(partition.node_to_district[i*4+j+1])
            #else
                print(partition.node_to_district[TH.parents_by_level[1][i*4+j+1]])
            end 
            println()
        end 
    end
end


# makes One PT swap attempt, non-parallel (i.e., there is only one processor, responsible for the entire PT situation)
function PT_swap!(
    TH::TemperingHierarchy, 
    partitions::Dict{Int, MultiLevelPartition}, 
    part_ID_to_temp::Dict{Int, Int}, 
    num_levels::Int, 
    pop_mins_maxs::Vector{Tuple{Float64, Float64}}, 
    measures, i, rng, betas,
    (no_merges, no_splits, almosts, succs)
)
    # (1). Set finers and coarsers 
    finers   = mod(i, 2) == 0 ? (1:2:(num_levels-1)) : (2:2:(num_levels-1)) 
    coarsers = mod(i, 2) == 0 ? (2:2:num_levels) : (3:2:num_levels)

    # (2). compute split node counts of each partition in a fnr level 
    split_nodes_by_level::Vector{Int} = [0 for _ = 1:(num_levels-1)]
    for part_ID in 1:num_levels 
        t = part_ID_to_temp[part_ID]
        if t in finers 
            split_nodes_by_level[t] = count_split_nodes(partitions[part_ID], TH.children_by_level[t+1])
        end
    end

    potential_swaps = Dict{Int, MultiLevelPartition}()

    # (3) compute the swap probabilities and stuff! 
    log_p1_factors = [0.0 for _ = 1:num_levels]
    log_p2_factors = [0.0 for _ = 1:num_levels]
    for part_ID in 1:num_levels 
        t = part_ID_to_temp[part_ID]
        if t in finers 
            fnr, csr = t, t+1
            f2c_partition = copy_partition(partitions[part_ID])
            log_p_up = swap_up(f2c_partition, TH, fnr, pop_mins_maxs, split_nodes_by_level[fnr], rng, betas[fnr])

            if log_p_up === nothing
                log_p1_factors[fnr] = -Inf 
                log_p2_factors[fnr] = -Inf 
                no_merges[fnr] += 1 
            else 
                cpartition_new = project_up_pctpair(f2c_partition, TH.graphs_by_level[csr], TH.parents_by_level[fnr]) 
                cur_energy = fully_compute_energy(partitions[part_ID], measures[fnr]) 
                new_energy = fully_compute_energy(cpartition_new, measures[csr])
                log_p1_factors[fnr] = log_p_up
                log_p2_factors[fnr] = (-new_energy + cur_energy)
                potential_swaps[part_ID] = cpartition_new
            end
        elseif t in coarsers
            fnr, csr = t-1, t
            c2f_partition = project_down_pctpair(partitions[part_ID], 
                TH.graphs_by_level[fnr], TH.children_by_level[csr])

            log_p_down = swap_down(c2f_partition, TH, fnr, pop_mins_maxs, split_nodes_by_level[fnr], rng, betas[fnr]) 

            if log_p_down === nothing 
                log_p1_factors[csr] = -Inf
                log_p2_factors[csr] = -Inf
                no_splits[csr] += 1
            else 
                cur_energy = fully_compute_energy(partitions[part_ID], measures[csr]) 
                new_energy = fully_compute_energy(c2f_partition, measures[fnr])
                log_p1_factors[csr] = log_p_down
                log_p2_factors[csr] = (-new_energy + cur_energy)
                potential_swaps[part_ID] = c2f_partition
            end
        else 
            # one of the ends, and we're not swapping it this time 
        end 
    end
    log_p_factors = log_p1_factors .+ log_p2_factors 

    # the finer indices will be true if the swap happens 
    swaps = [false for _=1:num_levels]
    for part_ID in 1:num_levels 
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

    for part_ID in 1:num_levels 
        t = part_ID_to_temp[part_ID]
        if t in finers && swaps[t]
            part_ID_to_temp[part_ID] = t+1
            partitions[part_ID] = potential_swaps[part_ID]
        elseif t in coarsers && swaps[t-1]
            part_ID_to_temp[part_ID] = t-1
            partitions[part_ID] = potential_swaps[part_ID]
        end
    end
end
