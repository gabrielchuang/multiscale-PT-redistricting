@enum Extensions REJECTION_COUNT ACCEPTANCE_COUNT TRANSITION_RATIOS TRANSITION_RATIOS2  

struct MultiLevelPartition
    num_dists::Int
    district_to_nodes::Vector{Mask}

    node_to_district::Dict{GlobalID, Int}
        # only has node_ID X as a key iff: 
        # 1) X is split between two districts, in which case value = -1, or 
        # 2) X is fully contained in district i, and X's parent is split. value = i
        # i.e., n2d is sparsely represented, with only the keys of the mask as keys in n2d. 
    graph::MultiLevelGraph
    subgraphs::Vector{MultiLevelGraph} # each subgraph's parent should be graph 
    dist_populations::Vector{Int}
    dist_areas::Vector{Real} 
    dist_perimeters::Vector{Real}

    adjacency::Array{Bool} # DxD, adjacency matrix for districts 

    # fix two districts, i and j. node X is a "shared node" of districts i and j iff: 
    # 0) i and j are adjacent 
    # 1) at least one descendant of X is in the mask for district i
    # 2) at least one descendant of X is in the mask for district j
    # 3) no descendants of X satisfy both 1) and 2). 
    # 4) the descendants in 1) and 2) are adjacent. 
    # so, e.g. if two districts each contain some of the 
    # precints of a county, that county is a shared node for those two districts. 
    # But if any precinct is split, then that precinct is a shared node, and the county is not. 

    shared_node::Array{GlobalID} 
    #2d DxD array, shared_node[i,j] is:
        # global_ID of shared node if i and j have exactly one shared node
        # 0 if i and j non-adjacent 
        # -1 if i and j share multiple nodes 
    # nb: shared_node[i,j] >0 is equiv to exists_hierarchical_tree(i,j). 

    adjacent_fine_neighbors::Array{Set{GlobalID}, 2}
    # 2D DxD array. NOT symmetric. 
    # each element [i,j] contains the precincts (or finest-level nodes, more generally), 
    # in district j that are adjacent to district i. 

    # a linking edge from d1 to d2 is an edge between sibling nodes, who belong to 
    # different districts, and whose parent is a shared_node. 
    linking_edge_attrs::Array{Union{Nothing, EdgeAttributes}}
    #2d DxD array, linking_edge_count[i,j] is the number of linking edges of dists i and j 
    # 0 if no linking edges or not exactly 1 shared node. 
    
    # not sure if needed for now 
    parent::Union{MultiLevelPartition,Nothing} # optional parent partition
    extensions::Dict{Union{Extensions, String},Any}
    proposed_extensions::Dict{Extensions,Any}
    favorite_pcts::Vector{Int64}

    articulation_points::Dict{Int, Set{GlobalID}}
end

""""""
function MultiLevelPartition(
    multi_level_graph::MultiLevelGraph,
    constraints::Dict,
    num_dists::Int;
    rng::AbstractRNG=PCG.PCGStateOneseq(UInt64),
    max_attempts::Int=100,
    only_cut_top_level_edge=false,
    first_comp_take=true
)#::Union{MultiLevelPartition, Nothing}
    compute_all_fine_neighbors(multi_level_graph) 
    for ii = 1:max_attempts
        partition, success = attemptConstructPartition(
            multi_level_graph, constraints, num_dists, rng, max_attempts, 
            only_cut_top_level_edge=only_cut_top_level_edge, first_comp_take=first_comp_take)
        if success

            append!(partition.favorite_pcts, [i*260 for i=1:10])

            partition.extensions[REJECTION_COUNT] = 
                Dict{String, Any}([
                    "rejection_prestep" => 0, 
                    "rejection_metropolis" => 0,
                    "rejection_noedge" => 0, 
                    "rejection_prestep_post" => 0,
                    "rejection_p0" => 0, 
                    "rejection_proposed_identical" => 0, 
                    "rejection_constraint" => 0, 
                    "continuous traversal constraint" => 0, 
                    "max split nodes by level" => 0, 
                    "max split nodes" => 0, 
                    "pop constraint" => 0, 
                    "intermediate rejections" => 0
                  ])
            partition.extensions[ACCEPTANCE_COUNT] = 0  
            partition.extensions[TRANSITION_RATIOS] = []
            partition.extensions[TRANSITION_RATIOS2] = []
            return partition
        end
    end
    println("failed max_attempts times")
    # failed max_attempts times; give up 
    return nothing 

    #throw(DomainError(success,
    #        "Failed to construct a random partition after "*
    #        string(max_attempts)*" attempts."))
end

function reconstruct_partition(
    graph::MultiLevelGraph, 
    leaf_nodes_by_district
)::MultiLevelPartition

end

function set_node_to_district!(
    node_to_district::Dict{GlobalID, Int}, 
    district_mask::Mask, 
    di::Int, 
    root::GlobalID 
)
    function set_district(node_ID) 
        if (district_mask[node_ID] === nothing || length(district_mask[node_ID]) == 0)
            node_to_district[node_ID] = di 
        else 
            node_to_district[node_ID] = -1 
        end 
    end 
    mask_preord_foreach(district_mask, root, set_district)
end 
function clear_node_to_district!(
    partition::MultiLevelPartition, 
    districts::Vector{Int}
)
    function clear_n2d(node_ID) 
       # println("deleting ", node_ID)
        delete!(partition.node_to_district, node_ID) 
    end 
    for d in districts 
        mask_preord_foreach(partition.district_to_nodes[d], partition.graph.root.global_ID, clear_n2d) 
    end 
end

# dependencies: d2n 
function set_node_to_district!(
    partition::MultiLevelPartition, 
    districts::Vector{Int}
)
    for d in districts
        set_node_to_district!(partition.node_to_district, partition.district_to_nodes[d], d, partition.graph.root.global_ID)
    end 
end 


""""""
function attemptConstructPartition(
    graph::MultiLevelGraph,
    constraints::Dict,
    num_dists::Int, 
    rng::AbstractRNG,
    max_attempts::Int;
    only_cut_top_level_edge::Bool = false, 
    first_comp_take = true # is this the first pop dev we're trying in expansion?
)
    max_attempts = 200
    district_to_nodes = Vector{Mask}(undef, num_dists)

    # remaining_nodes starts out as all the nodes 
    remaining_nodes::Mask = deepcopy(graph.mask) #Dict([graph.root.global_ID => nothing])  
    
    for district = 1:num_dists-1
        num_comp_fails = 0 

        found_district = false
        for attempt = 1:max_attempts
            subgraph = take_subgraph(graph, remaining_nodes)

            balance = (1, num_dists-district)
       
            multiscale_cuttable_tree = MultiScaleCuttableTree()
            edge_decs = Dict{GlobalID, Dict{GlobalID, Int}}()
            ext_decs = Dict{GlobalID, Dict{GlobalID, Int}}()

            #println("construct cuttable tree")
            construct_cuttable_tree!(multiscale_cuttable_tree, subgraph,
                constraints, subgraph.root.global_ID; 
                edge_decorators=edge_decs, ext_decorators=ext_decs, rng=rng, 
                balance=balance, 
                only_cut_top_level_edge=only_cut_top_level_edge)
            #flush(stdout)
            
            proposed_cut = cut_edge(multiscale_cuttable_tree, subgraph, rng,
                # if we only want to cut this "flatly", then we only want to cut edges who live 
                # under the root in the hierarchy. 
                restrict_node_cuts_to = only_cut_top_level_edge ? [graph.root.global_ID] : nothing)
            new_dist_masks, new_dist_pops, edge, prob_edge = proposed_cut

            #@show edge 

            if edge === nothing
                continue
            end 

            district_to_nodes[district] = new_dist_masks[1]

            picky_about_compactness = false 
            if picky_about_compactness && (first_comp_take || num_comp_fails < max_attempts - 10)
                d_finest_nodes = mask_all_finest_nodes(new_dist_masks[1])
                R = 6371
                area = sum([get_node_area(graph, x) for x in d_finest_nodes])
                min_lat = minimum([node_attributes(graph, x)["min_lat"] for x in d_finest_nodes])
                min_lon = minimum([node_attributes(graph, x)["min_lon"] for x in d_finest_nodes])
                max_lat = maximum([node_attributes(graph, x)["max_lat"] for x in d_finest_nodes])
                max_lon = maximum([node_attributes(graph, x)["max_lon"] for x in d_finest_nodes])
                box_area = abs(sin(pi/180 * max_lat) - sin(pi/180 * min_lat)) * abs(max_lon - min_lon) * pi/180 * R^2
                comp1 = area / box_area 

                d_finest_nodes = mask_all_finest_nodes(new_dist_masks[2])
                R = 6371
                area = sum([get_node_area(graph, x) for x in d_finest_nodes])
                min_lat = minimum([node_attributes(graph, x)["min_lat"] for x in d_finest_nodes])
                min_lon = minimum([node_attributes(graph, x)["min_lon"] for x in d_finest_nodes])
                max_lat = maximum([node_attributes(graph, x)["max_lat"] for x in d_finest_nodes])
                max_lon = maximum([node_attributes(graph, x)["max_lon"] for x in d_finest_nodes])
                box_area = abs(sin(pi/180 * max_lat) - sin(pi/180 * min_lat)) * abs(max_lon - min_lon) * pi/180 * R^2
                comp2 = area / box_area 

                threshold = 0.40
                if comp1 < threshold || comp2 < threshold
                    num_comp_fails += 1 
                    if num_comp_fails == max_attempts - 10 && !first_comp_take
                        println("hit max comp fails")
                    end
                    continue 
                end
            end
            
            if district == num_dists-1
                district_to_nodes[district+1] = new_dist_masks[2]
                tot_dists = num_dists
            else
                remaining_nodes_tmp = new_dist_masks[2]
                tot_dists = district
            end

            pop_const = constraints[PopulationConstraint]
            ideal_pop = (pop_const.min_pop + pop_const.max_pop)*0.5
            new_dist_pop_dev = round(100*(new_dist_pops[1]-ideal_pop)/ideal_pop, digits = 3)
            rem_dists = num_dists-district
            rem_avg_pop_dev = round(100*(new_dist_pops[2]/rem_dists-ideal_pop)/(ideal_pop), digits = 3)

            # check district constraints (not all)
            cur_districts = view(district_to_nodes, 1:tot_dists)
            # TODO(gtc) - ??? is this right 
            if !satisfies_constraint(ConstrainDiscontinuousTraversals(1), 
                take_subgraph(graph, new_dist_masks[1])) || 
                !satisfies_constraint(ConstrainDiscontinuousTraversals(1), 
                take_subgraph(graph, new_dist_masks[2]))
                continue

            elseif district < num_dists-1
                remaining_nodes = mask_intersection(remaining_nodes, remaining_nodes_tmp, graph)
                #print_mask(remaining_nodes, graph)
                if toggle_ENSURE()
                    remaining = take_subgraph(graph, remaining_nodes) 
                    #println("remaining pop", get_node_population(remaining, remaining.root.global_ID))
                    @ensure get_node_population(remaining, remaining.root.global_ID) == new_dist_pops[2]
                end 
            end

            # if success continue with remaining graph
            found_district = true
            #println("found district ", district )

            #print_mask(district_to_nodes[district], graph)
            if toggle_ENSURE() 
                district_i = take_subgraph(graph, district_to_nodes[district])
                #println("district i pop", get_node_population(district_i, district_i.root.global_ID))
                @ensure get_node_population(district_i, district_i.root.global_ID) == new_dist_pops[1]
            end 
            break
        end
        if !found_district
            return nothing, false
        end
    end 

    node_to_district = Dict{GlobalID, Int}() 
    dist_populations = Vector{Int}(undef, num_dists)
    dist_areas = Vector{Real}(undef, num_dists)
    dist_perimeters = Vector{Real}(undef, num_dists)

    subgraphs = Vector{MultiLevelGraph}(undef, num_dists)

    for di = 1:length(district_to_nodes) 
        # sets reverse map for each leaf in the mask 
        set_node_to_district!(node_to_district, district_to_nodes[di], di, graph.root.global_ID)

        district_mask::Mask = district_to_nodes[di]         
        district_subgraph = take_subgraph(graph, district_mask) 
        subgraphs[di] = district_subgraph
        dist_populations[di] = get_node_population(district_subgraph, district_subgraph.root.global_ID)
        dist_areas[di] = get_area(district_subgraph) 
    end 

    adjacency = falses((num_dists, num_dists))
    shared_node = zeros(GlobalID, (num_dists, num_dists))
    linking_edge_attrs = Array{Union{EdgeAttributes, Nothing}}(nothing, (num_dists, num_dists))

    parent = nothing
    extensions = Dict{String,Any}()
    proposed_extensions = Dict{String,Any}()
    adjacent_fine_neighbors = [Set() for _=1:num_dists, _=1:num_dists] 
    art_points = Dict{GlobalID, Set{GlobalID}}()

    partition = MultiLevelPartition(num_dists, district_to_nodes, node_to_district, graph,
        subgraphs, dist_populations, dist_areas, dist_perimeters, adjacency, shared_node, 
        adjacent_fine_neighbors, linking_edge_attrs, 
        parent, extensions, proposed_extensions, [], art_points)

    for di = 1:length(district_to_nodes) 
        dist_perimeters[di] = get_perimeter(partition, di)
        art_points[di] = articulation_points(partition, [di])
    end

    set_adjacency_and_shared_nodes!(partition)
    set_linking_edge_attrs!(partition)
    set_adjacent_fine_neighbors!(partition)
    return partition, true
end 

function partition_from_d2n(
    graph::MultiLevelGraph, 
    district_to_nodes, 
    num_dists; 
    extensions=Dict{String, Any}(),
    care_about_linking_edges = true
)
    node_to_district = Dict{GlobalID, Int}() 
    dist_populations = Vector{Int}(undef, num_dists)
    dist_areas = Vector{Real}(undef, num_dists)
    dist_perimeters = Vector{Real}(undef, num_dists)

    subgraphs = Vector{MultiLevelGraph}(undef, num_dists)

    for di = 1:num_dists  
        # sets reverse map for each leaf in the mask 
        set_node_to_district!(node_to_district, district_to_nodes[di], di, graph.root.global_ID)

        district_mask::Mask = district_to_nodes[di]         
        district_subgraph = take_subgraph(graph, district_mask) 
        subgraphs[di] = district_subgraph
        dist_populations[di] = get_node_population(district_subgraph, district_subgraph.root.global_ID)
        dist_areas[di] = get_area(district_subgraph) 
    end 

    adjacency = falses((num_dists, num_dists))
    shared_node = zeros(GlobalID, (num_dists, num_dists))
    linking_edge_attrs = Array{Union{EdgeAttributes, Nothing}}(nothing, (num_dists, num_dists))

    parent = nothing
    proposed_extensions = Dict{String,Any}()
    adjacent_fine_neighbors = [Set() for _=1:num_dists, _=1:num_dists] 
    art_points = Dict{GlobalID, Set{GlobalID}}()

    partition = MultiLevelPartition(num_dists, district_to_nodes, node_to_district, graph,
        subgraphs, dist_populations, dist_areas, dist_perimeters, adjacency, shared_node, 
        adjacent_fine_neighbors, linking_edge_attrs, 
        parent, extensions, proposed_extensions, [], art_points)

    for di = 1:length(district_to_nodes) 
        dist_perimeters[di] = get_perimeter(partition, di)
        art_points[di] = articulation_points(partition, [di])
    end

    set_adjacency_and_shared_nodes!(partition)
    if care_about_linking_edges
        set_linking_edge_attrs!(partition)
    end
    set_adjacent_fine_neighbors!(partition)
    return partition
end

# pretty expensive, but we should only do this once (at creation of partition). 
function set_adjacent_fine_neighbors!(partition::MultiLevelPartition) 
    for e in finest_level_edges(partition.graph)
        if !is_in_graph(partition.graph, e.id1) || !is_in_graph(partition.graph, e.id2) 
            continue 
        end
        d1 = look_up_district(partition, e.id1)
        d2 = look_up_district(partition, e.id2)
        if d1 != d2 
            push!(partition.adjacent_fine_neighbors[d1,d2], e.id2)
            push!(partition.adjacent_fine_neighbors[d2,d1], e.id1)
        end
    end
end


# iterating throught the multilevelneighbors of curr_node to find adjacencies and shared_nodes, 
# setting them in adjacency and shared_node respectively. 
function add_nbr_districts_and_shared_nodes!(
    adjacency::Matrix{Bool}, 
    shared_node::Union{Nothing, Matrix{Int}},
    partition::MultiLevelPartition, 
    curr_dist::Int, 
    curr_nbr::MultiLevelNeighbor, 
    curr_node::GlobalID
) 
    nbr_dist = look_up_district(partition, curr_nbr.global_ID)
    if nbr_dist == -1 
        for finer_nbr in values(curr_nbr.finer_neighbors) 
            add_nbr_districts_and_shared_nodes!(adjacency, shared_node, partition, curr_dist, finer_nbr, curr_node)
        end 
    elseif nbr_dist !== curr_dist
        adjacency[nbr_dist, curr_dist] = true 
        adjacency[curr_dist, nbr_dist] = true 

        if shared_node !== nothing 
            potentially_shared_node = least_common_ancestor(partition.graph, curr_node, curr_nbr.global_ID)

            if shared_node[nbr_dist, curr_dist] == potentially_shared_node
                #we already set the shared node for this pair of districts, and it's the same  
            elseif shared_node[nbr_dist, curr_dist] == 0 || #we haven't set the shared node yet, or 
                is_descendant(partition.graph, potentially_shared_node, shared_node[nbr_dist, curr_dist])
                # we set it already to an ancestor of this one, so we can update it to be finer 
                shared_node[nbr_dist, curr_dist] = potentially_shared_node
                shared_node[curr_dist, nbr_dist] = potentially_shared_node
            elseif shared_node[nbr_dist, curr_dist] > 0 && # we already set it, and 
                !is_descendant(partition.graph, shared_node[nbr_dist, curr_dist], potentially_shared_node)
                # the thing we set it to isn't a descendant of this possibility, so there are multiple 
                # "shared nodes", which means that this pair of dists aren't eligible to be combined. 
                shared_node[nbr_dist, curr_dist] = -1
                shared_node[curr_dist, nbr_dist] = -1
            end 
        end 
    end 
end 

function clear_adjacency!(partition::MultiLevelPartition, districts::Vector{Int})
    for d in districts 
        partition.adjacency[d,:] .= 0 
        partition.adjacency[:,d] .= 0 
    end 
    @assert issymmetric(partition.adjacency)
end 

# dependencies: n2d, fine_neighbors
function set_adjacency_and_shared_nodes!(partition::MultiLevelPartition, changed_nodes=1:partition.num_dists)
    n2d = partition.node_to_district

    for di in changed_nodes 
        district_mask = partition.district_to_nodes[di]
        function add_adjacencies(node_ID)
            @ensure haskey(n2d, node_ID) && n2d[node_ID] != -1 
            for multilevel_nbr in fine_neighbors(partition.graph, node_ID)
                add_nbr_districts_and_shared_nodes!(partition.adjacency, partition.shared_node, 
                    partition, n2d[node_ID], multilevel_nbr, node_ID)
            end 
        end 

        mask_preord_foreach_leaf_only(district_mask, partition.graph.root.global_ID, add_adjacencies)
    end
    @ensure issymmetric(partition.adjacency)
end 

function clear_shared_nodes!(partition::MultiLevelPartition, districts::Vector{Int}) 
    for d in districts 
        partition.shared_node[d,:] .= 0 
        partition.shared_node[:,d] .= 0 
    end 
end 

function clear_linking_edge_attrs!(partition::MultiLevelPartition, districts::Vector{Int}) 
    for d in districts 
        partition.linking_edge_attrs[d,:] .= nothing
        partition.linking_edge_attrs[:,d] .= nothing
    end 
end 

function count_finest_level_cut_edges(partition::MultiLevelPartition) 
    count = 0 
    for edge in finest_level_edges(partition.graph)
        v1, v2 = edge.id1, edge.id2
        @ensure look_up_district(partition, v1) != -1 
        @ensure look_up_district(partition, v2) != -1 

        if look_up_district(partition, v1) != look_up_district(partition, v2) 
            count += 1 
        end 
    end
    return count 
end

function count_finest_level_cut_edges(partition, _) 
    count_finest_level_cut_edges(partition)
end 


# returns either the total weight of linking edges (if just_count=true) 
# or (total weight, vector of UndirectedEdge{GlobalID}) 
function collect_all_linking_edges(partition, d1::Int, d2::Int; just_count=false) 
    @ensure partition.shared_node[d1,d2] > 0 
    g, vmap = child_graph(partition.graph, partition.shared_node[d1,d2])
    linking_edges = just_count ? nothing : [] 
    total_attrs = EdgeAttributes(0,0)

    #initial_possibilities = [(vmap[src(e)], vmap[dst(e)]) for e in edges(g)]

    function accumulate_linking_edges(potential_edges, pre_vmap=false)
        for e in potential_edges 
            e1, e2 = -1, -1
            if pre_vmap 
                e1, e2 = vmap[src(e)], vmap[dst(e)]
            else
                e1, e2 = e 
            end
            associated_districts = (look_up_district(partition, e1), look_up_district(partition, e2))
            if associated_districts == (d1, d2) || associated_districts == (d2, d1) 
                if !just_count
                    push!(linking_edges, UndirectedEdge(e1, e2)) 
                end 
                total_attrs = combine_edge_attrs(total_attrs, get_edge_attributes(partition.graph, e1, e2))
            elseif associated_districts == (d1, -1) || associated_districts == (d2, -1) 
                fine_neighbors(partition.graph, e1)
                potential_e2s = keys(partition.graph.all_neighbors_map[e1][e2].finer_neighbors)
                new_potential_edges = [(e1, x) for x in potential_e2s]
                accumulate_linking_edges(new_potential_edges)
            elseif associated_districts == (-1, d1) || associated_districts == (-1, d2) || 
                    associated_districts == (-1, -1)
                fine_neighbors(partition.graph, e2)
                potential_e1s = keys(partition.graph.all_neighbors_map[e2][e1].finer_neighbors)
                new_potential_edges = [(e2, x) for x in potential_e1s]
                accumulate_linking_edges(new_potential_edges)
            end
        end
    end

    accumulate_linking_edges(edges(g), true)

    if just_count 
        return total_attrs 
    else 
        return (total_attrs, linking_edges) 
    end 
end

function set_linking_edge_attrs!(partition::MultiLevelPartition, changed_nodes=1:partition.num_dists) 
    for i in changed_nodes  
        for j=1:partition.num_dists 
            if partition.shared_node[i,j] > 0 
                linking_edge_attrs = collect_all_linking_edges(partition, i, j, just_count=true) 
                @assert linking_edge_attrs.weight != 0 
                partition.linking_edge_attrs[i,j] = linking_edge_attrs
                partition.linking_edge_attrs[j,i] = linking_edge_attrs
            end 
        end 
    end 
end 

function log_hierarchical_forest_count(P::MultiLevelPartition) 
    for s in P.subgraphs
        precompute_log_tree_counts!(s) 
    end 
    return sum([log_hierarchical_tree_count(s) for s in P.subgraphs])
end 

function look_up_district(partition::MultiLevelPartition, node_ID::GlobalID) 
    if haskey(partition.node_to_district, node_ID) 
        return partition.node_to_district[node_ID] 
    end
    return look_up_district(partition, partition.graph.node_parents[node_ID])
end 

function write_partition_to_file(partition::MultiLevelPartition, output_file::String, dists=1:partition.num_dists)
    out = []
    G = partition.graph

    num_named_levels = length(node_name_exclude_quotient_nodes(
        partition.graph, partition.graph.ID_to_node[1].global_ID))

    for dist in dists
        function f(G, k)
            node_name = node_name_exclude_quotient_nodes(G, k.global_ID)
            if length(node_name) == num_named_levels
                s = join(node_name, ",")
                push!(out, "\""*s*"\",\""*string(dist)*"\"")
            end
        end
    
        for (k,v) in partition.district_to_nodes[dist] 
            if v === nothing 
                if G.ID_to_node[k].node_name !== nothing 
                    node_name = node_name_exclude_quotient_nodes(G, k)
                    s = join(node_name, ",")
                    push!(out, "\""*s*"\",\""*string(dist)*"\"")
                elseif length(node_name_exclude_quotient_nodes(G, k)) < num_named_levels
                    preord_foreach(partition.subgraphs[dist], k, f)
                end
            end
        end
    end
    s = "node_name,district\n"*join(out, "\n")
    
    open(output_file, "w") do file 
        write(file, s)
    end 
end 

function write_partition_to_file_flat(partition::MultiLevelPartition, output_file::String) 
    out = []
    for node in finest_level_nodes(partition.graph) 
        d = look_up_district(partition, node)
        push!(out, string(node) * ", " * string(d))
    end
    s = "node_name,district\n"*join(out, "\n")
    
    open(output_file, "w") do file 
        write(file, s)
    end 
end

