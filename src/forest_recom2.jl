#using GraphPlot, Compose, ColorTypes 
#import Cairo, Fontconfig 

# returns a uniform choice out of the pairs of districts who are 
# 1) adjacent, and 
# 2) can have a hierarchical tree drawn on them (i.e. shared_node > 0), and 
# 3) have both districts in districts_to_choose_from
function sample_combinable_districts_randomly(
    partition::MultiLevelPartition,
    rng::AbstractRNG, 
    districts_to_choose_from::Vector{Int} 
)
    adj_dists = Vector{Pair{Int,Int}}(undef, 0)
    for di = 1:length(districts_to_choose_from)-1
        for dj = di+1:length(districts_to_choose_from)
            if partition.shared_node[districts_to_choose_from[di], districts_to_choose_from[dj]] > 0
                push!(adj_dists, Pair(districts_to_choose_from[di], districts_to_choose_from[dj]))
            end
        end
    end
    dist_pair = rand(rng, adj_dists)
    return dist_pair, 1.0/length(adj_dists)
end 

function sample_any_districts_randomly(
    partition::MultiLevelPartition,
    rng::AbstractRNG
)
    adj_dists = Vector{Pair{Int,Int}}(undef, 0)
    for di = 1:partition.num_dists-1
        for dj = di+1:partition.num_dists
            if partition.adjacency[di, dj]
                push!(adj_dists, Pair(di, dj))
            end
        end
    end
    dist_pair = rand(rng, adj_dists)
    return dist_pair, 1.0/length(adj_dists)
end 

# pick two random districts. return the district numbers, 1/num adjacencies, and the 
# induced subgraph from the union of these two dists. 
function sample_subgraph(
    partition::MultiLevelPartition,
    rng::AbstractRNG, 
    distrits_to_choose_from::Vector{Int})
    proposed_dists = sample_combinable_districts_randomly(partition, rng, distrits_to_choose_from)
    #proposed_dists = sample_any_districts_randomly(partition, rng)

    (di, dj), prob_of_dists = proposed_dists

    merged_mask = mask_union(partition.district_to_nodes[di], 
        partition.district_to_nodes[dj], 
        partition.graph)

    subgraph = take_subgraph(partition.graph, merged_mask)
    return [di, dj], prob_of_dists, subgraph
end

# pick a linking edge between d1 and d2 (weighted by edge weight) 
function sample_linking_edge(
    partition::MultiLevelPartition,
    d1::Int, d2::Int, rng::AbstractRNG
)::UndirectedEdge{GlobalID}
    @ensure partition.adjacency[d1, d2] 
    @ensure partition.shared_node[d1, d2] > 0 
    @ensure partition.linking_edge_attrs[d1,d2] !== nothing 

    shared_node_child_graph, vmap = child_graph(partition.graph, partition.shared_node[d1, d2])
    choice = rand(rng) * partition.linking_edge_attrs[d1, d2].weight 
    cumsum = 0

    (_, linking_edges) = collect_all_linking_edges(partition, d1, d2) 
    for e in linking_edges 
        cumsum += get_edge_weight(partition.graph, e.id1, e.id2) 
        if cumsum >= choice 
            return e
        end
    end 

    write_partition_to_file(partition, "viz/err.csv", [d1,d2])
    @assert false 
end

#refine_linking_edge() takes an edge (a,b) and returns an edge (c,d) such that: 
# - c and d are both leaves of the graph (i.e. cblocks) 
# - c is a child of a, and b is a child of b (or vice versa) 
# - c and d are joined by an edge 
function refine_linking_edge(
    G::MultiLevelGraph,
    node1::GlobalID, 
    node2::GlobalID,
    rng::AbstractRNG
)
    if length(get_children_IDs(G, node1)) == 0 
        if length(get_children_IDs(G, node2)) == 0 
            return UndirectedEdge(node1, node2) 
        else 
            return refine_linking_edge(G, node2, node1, rng) 
        end 
    end 
    #@assert node1 has children 

    node1_refinement = refine_node(G, node1, node2, rng) 
    return refine_linking_edge(G, node1_refinement, node2, rng)
end 


function sample_merged_tree(
    changed_districts::Vector{Int},
    partition::MultiLevelPartition,
    constraints::Dict,
    rng::AbstractRNG
)
    # prob_selecting_old_edge_to_remove
    D1, D2 = changed_districts
    edge = sample_linking_edge(partition, D1, D2, rng)
    @assert edge !== nothing 

    subgraph1 = partition.subgraphs[D1]
    subgraph2 = partition.subgraphs[D2]

    specified_edge = refine_linking_edge(partition.graph, edge.id1, edge.id2, rng)

    multiscale_cuttable_tree1 = MultiScaleCuttableTree()
    multiscale_cuttable_tree2 = MultiScaleCuttableTree()
    pop1 = partition.dist_populations[D1]
    pop2 = partition.dist_populations[D2]

    edge_end_1 = is_in_graph(subgraph1, specified_edge.id1) ? specified_edge.id1 : specified_edge.id2
    edge_end_2 = is_in_graph(subgraph1, specified_edge.id1) ? specified_edge.id2 : specified_edge.id1

    exterior_decorators1 = Dict{GlobalID,Dict{GlobalID,Int}}(
        edge_end_1 => Dict(edge_end_2 => pop2))
    exterior_decorators2 = Dict{GlobalID,Dict{GlobalID,Int}}(
        edge_end_2 => Dict(edge_end_1 => pop1))

    construct_cuttable_tree!(multiscale_cuttable_tree1, subgraph1, constraints,
        ext_decorators=exterior_decorators1, rng=rng)
    construct_cuttable_tree!(multiscale_cuttable_tree2, subgraph2, constraints,
        ext_decorators=exterior_decorators2, rng=rng)

    total_cuttable_edges = (count_cuttable_edges(multiscale_cuttable_tree1)
        + count_cuttable_edges(multiscale_cuttable_tree2) 
        + 1) #because the actual edge...? 
    edge_prob = 1.0 / total_cuttable_edges
    @assert total_cuttable_edges > 0 


    old_tree = (multiscale_cuttable_tree1, multiscale_cuttable_tree2, specified_edge)
    return old_tree, edge_prob
end

# how many adjacencies are there to `district_mask` in `partition`? 
function count_adjacent_dists(
    partition::MultiLevelPartition, 
    district::Int, 
    district_mask::Mask
)::Int 
    adj_dists = Matrix{Bool}(falses((partition.num_dists, partition.num_dists)))

    function add_adjacencies(node_ID)
        for multilevel_nbr in fine_neighbors(partition.graph, node_ID)
            add_nbr_districts_and_shared_nodes!(adj_dists, nothing, partition, 
                district, multilevel_nbr, node_ID)
        end 
    end 
    mask_preord_foreach_leaf_only(district_mask, partition.graph.root.global_ID, add_adjacencies)
    @ensure !adj_dists[district, district] 
    @ensure issymmetric(adj_dists) 
    return count(adj_dists)/2  
end 

# suppose we had `partition`, but with d1_mask and d2_mask replacing d1 and d2. 
# when we pick two random adjacent districts, what's the probability that we pick d1 and d2? 
function get_prob_of_new_dists(
    partition::MultiLevelPartition,
    d1::Int, d2::Int, 
    d1_mask::Mask, d2_mask::Mask, 
)::Real 
    @ensure partition.adjacency[d1,d2]
    num_pairs_between_non_d1_d2_districts = 
        count(partition.adjacency) / 2 - #adjacency in the OG partition 
        count(partition.adjacency[d1,:]) - count(partition.adjacency[d2,:]) + 1
        # +1 since the adjacency from d1 to d2 is subtracted twice 
    num_adjacencies_for_d1 = count_adjacent_dists(partition, d1, d1_mask) 
    num_adjacencies_for_d2 = count_adjacent_dists(partition, d2, d2_mask) 
    num_pairs_involving_d1_d2 = num_adjacencies_for_d1 + num_adjacencies_for_d2 - 1 
    # the adjacency from d1 to d2 is counted twice 

    total_pairs = num_pairs_between_non_d1_d2_districts + num_pairs_involving_d1_d2
    return 1.0 / total_pairs 
end

function log_linking_edge_prob_tree_space(
    partition::MultiLevelPartition,
    d1::Int, d2::Int
)
    ans = - (
        sum(log.([x.weight for x in filter(x -> x !== nothing, partition.linking_edge_attrs[d1,:])]))
        + sum(log.([x.weight for x in filter(x -> x !== nothing, partition.linking_edge_attrs[d2,:])]))
        - log(partition.linking_edge_attrs[d1,d2].weight))
    return ans 
    # the weight from d1 to d2 is counted twice, so it's subtracted out.  
end

function log_linking_edge_prob_tree_space(
    partition::MultiLevelPartition, 
)
    ans = 0 
    for cd in 1:partition.num_dists 
        ans += sum(log.([x.weight for x in filter(x -> x !== nothing, partition.linking_edge_attrs[cd,:])]))
    end
    return -ans / 2.0 #every edge is counted twice 
end

function take_snapshot(
    partition::MultiLevelPartition, 
    changed_dists::Vector{Int}
)
    old_masks = [partition.district_to_nodes[i] for i in changed_dists]
    old_pops = [partition.dist_populations[i] for i in changed_dists]    
    old_perims = [partition.dist_perimeters[i] for i in changed_dists]    
    old_areas = [partition.dist_areas[i] for i in changed_dists]    
    old_adjacency = copy(partition.adjacency) 
    old_shared_nodes = copy(partition.shared_node) 
    old_linking_edge_attrs = copy(partition.linking_edge_attrs)

    return old_masks, old_pops, old_perims, old_areas, old_adjacency, old_shared_nodes, old_linking_edge_attrs
end 

function take_snapshot_and_update!(
    partition::MultiLevelPartition, 
    changed_dists::Vector{Int}, 
    changed_masks::Vector{Mask}, 
    changed_pops::Vector{Int}
) 

    snapshot = take_snapshot(partition, changed_dists)

    clear_node_to_district!(partition, changed_dists)
    clear_adjacency!(partition, changed_dists)
    clear_shared_nodes!(partition, changed_dists)
    clear_linking_edge_attrs!(partition, changed_dists)

    for (i, di) in enumerate(changed_dists)
        partition.district_to_nodes[di] = changed_masks[i]
        partition.subgraphs[di] = take_subgraph(partition.graph, changed_masks[i])
        partition.dist_populations[di] = changed_pops[i] 
    end
    set_node_to_district!(partition, changed_dists)
    set_adjacency_and_shared_nodes!(partition, changed_dists)
    set_linking_edge_attrs!(partition, changed_dists)

    for (ii,cd) in enumerate(changed_dists)
        partition.dist_perimeters[cd] = get_perimeter(partition, cd)
        partition.dist_areas[cd] = get_area(partition.subgraphs[cd])
    end

    return changed_dists, snapshot...
end 

function restore_from_snapshot!(
    partition::MultiLevelPartition, 
    old_dists, old_masks, old_pops, old_perims, old_areas, old_adjacency, old_shared_nodes, old_linking_edge_attrs
)
    clear_node_to_district!(partition, old_dists)
    clear_adjacency!(partition, old_dists)
    clear_shared_nodes!(partition, old_dists)
    clear_linking_edge_attrs!(partition, old_dists)

    for (i, od) in enumerate(old_dists) 
        partition.district_to_nodes[od] = old_masks[i]
        partition.dist_populations[od] = old_pops[i] 
        partition.dist_perimeters[od] = old_perims[i] 
        partition.dist_areas[od] = old_areas[i] 
        partition.subgraphs[od] = take_subgraph(partition.graph, old_masks[i])
    end 
    set_node_to_district!(partition, old_dists)
    copyto!(partition.adjacency, old_adjacency)
    copyto!(partition.shared_node, old_shared_nodes)
    copyto!(partition.linking_edge_attrs, old_linking_edge_attrs)
end 



function check_prestep_constraints(
    partition::MultiLevelPartition, 
    constraints, 
    subgraph::MultiLevelGraph, 
    changed_districts
)
    if partition.shared_node[changed_districts...] == 0
        @assert false 
        return false 
    elseif haskey(constraints, ConstrainDiscontinuousTraversals) && 
        !satisfies_constraint(constraints[ConstrainDiscontinuousTraversals], subgraph)
        #write_partition_to_file(partition, "viz/failed_prestep.csv")
        return false 
    end
    return true 
end 

""""""
function forest_recom2!(
    partition::MultiLevelPartition,
    measure::Measure,
    constraints::Dict,
    rng::AbstractRNG, 
    i::Int;
    districts_to_choose_from=collect(1:partition.num_dists)
)

    changed_districts, prob_of_dists, subgraph = sample_subgraph(partition, rng, districts_to_choose_from)

    d1, d2 = changed_districts 
    @assert partition.adjacency[d1,d2]

    if !check_prestep_constraints(partition, constraints, subgraph, changed_districts)
        if haskey(partition.extensions, REJECTION_COUNT)
            partition.extensions[REJECTION_COUNT]["rejection_prestep"]+=1
        end
        return 0, nothing
    end 

    multiscale_cuttable_tree = MultiScaleCuttableTree()
    construct_cuttable_tree!(multiscale_cuttable_tree, subgraph, constraints,
        subgraph.root.global_ID; rng=rng)

    proposed_cut = cut_edge(multiscale_cuttable_tree, subgraph, rng)
    new_dist_masks, new_dist_pops, edge, prob_edge = proposed_cut

    if edge === nothing
        partition.extensions[REJECTION_COUNT]["rejection_noedge"]+=1
        #return 0, nothing
    end    

    num_fails = 0 
    while edge === nothing && num_fails < 100
        num_fails += 1
        changed_districts, prob_of_dists, subgraph = sample_subgraph(partition, rng, districts_to_choose_from)

        d1, d2 = changed_districts 


        if !check_prestep_constraints(partition, constraints, subgraph, changed_districts)
            if haskey(partition.extensions, REJECTION_COUNT)
                partition.extensions[REJECTION_COUNT]["rejection_prestep"]+=1
            end
            return 0, nothing
        end 

        multiscale_cuttable_tree = MultiScaleCuttableTree()
        construct_cuttable_tree!(multiscale_cuttable_tree, subgraph, constraints,
            subgraph.root.global_ID; rng=rng)

        proposed_cut = cut_edge(multiscale_cuttable_tree, subgraph, rng)
        new_dist_masks, new_dist_pops, edge, prob_edge = proposed_cut

        if edge === nothing
            partition.extensions[REJECTION_COUNT]["rejection_noedge"]+=1
            #return 0, nothing
        end
    end 

    if num_fails == 100
        println("failed too many times")
        return 0, nothing 
    end 
    
    if !satisfies_constraints(partition, constraints, (changed_districts, new_dist_masks, new_dist_pops))
        partition.extensions[REJECTION_COUNT]["rejection_constraint"] += 1
        return 0, nothing
    end
    
    d1_mask, d2_mask = new_dist_masks 
    d1_pop, d2_pop = new_dist_pops

    if d1_mask == partition.district_to_nodes[d1] || d1_mask == partition.district_to_nodes[d2] 
        partition.extensions[REJECTION_COUNT]["rejection_proposed_identical"]+=1
        return (0, nothing) 
    end

    _, prob_old_edge = sample_merged_tree(changed_districts, 
        partition, constraints, rng)

    if prob_old_edge == 0
        println("couldn't draw old tree -- problem")
        @assert false
    end

    if measure.gamma == 0
        snapshot = take_snapshot_and_update!(partition, changed_districts, new_dist_masks, new_dist_pops) 

        # this is the BACKWARD prob: i.e., the probability that, given the new partition, 
        # d1 and d2 were chosen to be merge/splitted. 

        _, prob_of_new_dists = sample_combinable_districts_randomly(partition, rng, districts_to_choose_from)
        partialFwdPrpProb = prob_of_dists*prob_edge
        partialBckPrpProb = prob_of_new_dists*prob_old_edge
        p = partialBckPrpProb/partialFwdPrpProb

        restore_from_snapshot!(partition, snapshot...)
 
    elseif measure.gamma != 0 #&& p != 0 
        # nodes that are split between d1 and d2 in the proposal *or* the current partition
        # these are in fact just the ancestors of the shared node in each partition 
        shared_node_ancestors::Vector{GlobalID} = [] 
        # the ancestors of the shared node in the current partition 
        add_ancestors!(shared_node_ancestors, subgraph, partition.shared_node[d1, d2])

        # compute the log linking edge/tree count ratios 
        log_l_edge_prob_cur = log_linking_edge_prob_tree_space(partition, d1, d2)
       
        snapshot = take_snapshot_and_update!(partition, changed_districts, new_dist_masks, new_dist_pops) 
        # the ancestors of the shared node in the proposal partition 
        add_ancestors!(shared_node_ancestors, subgraph, partition.shared_node[d1, d2])
        shared_node_ancestors = unique!(shared_node_ancestors)

        _, prob_of_new_dists = sample_combinable_districts_randomly(partition, rng, districts_to_choose_from)
        partialFwdPrpProb = prob_of_dists*prob_edge
        partialBckPrpProb = prob_of_new_dists*prob_old_edge
        p = partialBckPrpProb/partialFwdPrpProb
        
        log_l_edge_prob_proposed = log_linking_edge_prob_tree_space(partition, d1, d2)
        log_tree_count_proposed = get_log_tree_counts(partition, subgraph, d1, d2, shared_node_ancestors)
        
        if haskey(constraints, ConstrainDiscontinuousTraversals) && 
            (!satisfies_constraint(constraints[ConstrainDiscontinuousTraversals], partition.subgraphs[d1]) || 
                !satisfies_constraint(constraints[ConstrainDiscontinuousTraversals], partition.subgraphs[d2]))
            log_tree_count_proposed = -Inf
        end 

        restore_from_snapshot!(partition, snapshot...)

        log_tree_count_cur = get_log_tree_counts(partition, subgraph, d1, d2, shared_node_ancestors)

        log_linking_edge_ratio = log_l_edge_prob_proposed - log_l_edge_prob_cur
        log_tree_count_ratio = log_tree_count_proposed - log_tree_count_cur 

        multiplier = exp(measure.gamma*(log_linking_edge_ratio + log_tree_count_ratio))

        p *= multiplier
    end

    return p, ([d1, d2], new_dist_masks, new_dist_pops, edge)
end

function build_forest_recom2(
    constraints::Dict
)
    f(p, m, r, i) = forest_recom2!(p, m, constraints, r, i)
    return f
end
