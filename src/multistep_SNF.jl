# the proposal for multistep-SNF. 
# given a starting state, makes a bunch of tempered SNFs in a given direction 
# e.g., a given cycle, then returns the accept prob 
# (based on the forward/reverse transition probs and the measure at each end)

function multistep_SNF(
    partition::MultiLevelPartition, 
    measure::Measure, 
    constraints::Dict, 
    rng::AbstractRNG, 
    num_steps::Int = 100 #eventually replace this with NUTS 
) 

    # pick a district cycle (or dist pair) to evolve and generate the respective filters 
    # todo: how to do this? do we still want to pick using the county-split tree on the 
    #   district graph? having >= n-1 split counties to use among a n-district rotate 
    #   is probably desirable
    #   it might make sense to just maintain a set of split counties (in the partition)
    #cycle, p_cycle = sample_district_cycle(partition, rng) #e.g., [dist 1 -> dist 2, 2 -> 3, 3 -> 1]
    tri, p_tri = pick_random_triangle(partition, rng)
    if rand(rng, [1,2]) == 1 
        tri = reverse(tri)
    end 
    tri = collect(tri)
    #tri = [12, 14, 9]
    #tri = [4,9,10]
    rev_tri = reverse(tri)

    filter::Vector{Tuple{Int64, Int64}} = make_filter_from_cycle(collect(tri))
    rev_filter::Vector{Tuple{Int64, Int64}} = reverse([(y,x) for (x,y) in filter])

    @show filter 

    forward_probs = [] 
    backward_probs = [] 

    forward_npcs = [] 
    backward_npcs = []

    sec_forward_probs = []
    sec_backward_probs = []

    fwdranks = [] 
    bwdranks = []

    flips = []
    merged_areas = lie_about_articulation_points(partition, tri, rng)
    @show merged_areas

    # compute the measure on the initial state 
    #p_initial = exp(-fully_compute_energy(partition, measure))
    track = Dict([
        "pop" => ([], [], []), 
        "comp" => ([], [], []), 
        "ctysplits" => [], 
        "energy" => []
    ])

    for i=1:num_steps 
        #@show i 
        # pick a random filter? (e.g., a random dist pair in the cycle)
        #which_filter::Vector{Tuple{Int64, Int64}} = [rand(rng, filter)] 
        
        # make a SNF with this filter and compute fwd/rev probabilities 
        prob, reject_update, accept_update, pct, dist_to, npcs_fwd, sec_prob, which_flip, fwdrank, which_flip2 = 
            single_node_flip(partition, measure, constraints, rng, filter; 
            merged_nodes=merged_areas, i=i)
        dist_from = look_up_district(partition, pct)
        accept_update(nothing) 

        push!(forward_probs, log(prob)) 
        push!(forward_npcs, npcs_fwd)
        push!(sec_forward_probs, log(sec_prob))
        push!(flips, which_flip)
        
        if i > 1000000
            @printf("move: %s, p_fwd = %.3f, pct = %s, dist_from = %d, dist_to = %d\n", 
                node_name(partition.graph, pct), prob, pct, dist_from, dist_to);
        end 
        #rev_which::Vector{Tuple{Int64, Int64}} = [(which_filter[1][2], which_filter[1][1])]
        backward_prob, npcs_bwd, sec_backward_prob, bwdrank = 
            single_node_flip(partition, measure, constraints, rng, 
            rev_filter, just_prob_of_flip=(dist_from, pct, !which_flip2); merged_nodes=merged_areas, i=i)

        @assert backward_prob != 0 

        push!(backward_probs, log(backward_prob)) 
        push!(backward_npcs, npcs_bwd)
        push!(sec_backward_probs, log(sec_backward_prob))

        push!(fwdranks, fwdrank)
        push!(bwdranks, bwdrank) 

        push!(track["pop"][1], partition.dist_populations[tri[1]])
        push!(track["pop"][2], partition.dist_populations[tri[2]])
        push!(track["pop"][3], partition.dist_populations[tri[3]])
        push!(track["comp"][1], partition.dist_perimeters[tri[1]]^2/partition.dist_areas[tri[1]])
        push!(track["comp"][2], partition.dist_perimeters[tri[2]]^2/partition.dist_areas[tri[2]])
        push!(track["comp"][3], partition.dist_perimeters[tri[3]]^2/partition.dist_areas[tri[3]])
        push!(track["ctysplits"], count_split_counties_multiplicity(partition))
        push!(track["energy"], fully_compute_energy(partition, measure))        

        if mod(i, 100) == 0 || i == num_steps-1
            write_partition_to_file(partition, "viz/step"*string(i)*".csv", tri)
        end
    end 
    # compute the measure on the final state x
    #p_final = exp(-fully_compute_energy(partition, measure))

    return forward_probs, backward_probs, track, forward_npcs, backward_npcs, 
        sec_forward_probs, sec_backward_probs, flips, fwdranks, bwdranks

    total_log_fwd_prob = sum(forward_probs) 
    total_log_bck_prob = sum(backward_probs) 

    # compute the acceptance prob 
    A = min(1, p_final / p_initial * exp(total_log_fwd_prob - total_log_bck_prob))

    # another question: how to revert? (i.e., what to do if we decide, *outside* this 
    # function, to reject the update?) 
    # the way I did it before was to return two functions 
    # (kinda like how s_n_f2 does now) for accept and reject; reject would maybe hold 
    # a deep copy of the initial state, or a subset of n2d or d2n (former reqs more space, 
    # latter reqs more computation?)
    return A, []

end

function nodes_bounded_by(
    partition::MultiLevelPartition, 
    rng;
    dists::Vector{Int}=1:partition.num_dists, 
    bounding_pts::Set{GlobalID}=[], 
)::Vector{Tuple{Set{GlobalID}, Set{GlobalID}}}
    res = []
    vxs = Set{GlobalID}()
    for i in dists 
        preord_foreach(partition.subgraphs[i], partition.graph.root.global_ID, 
            (G,x) -> length(get_children_IDs(G, x)) == 0 ? push!(vxs, x.global_ID) : nothing)
    end 
    while length(vxs) > 0 
        vx = rand(rng, vxs) 
        visited = Set() 
        this_area_bounding_pts = Set() 
        set_nodes_bounded_by!(partition, dists, bounding_pts, visited, this_area_bounding_pts, vx)
        push!(res, (visited, this_area_bounding_pts))
        setdiff!(vxs, visited) 
    end
    return res 
end 

# the set of visited nodes INCLUDES bounding points. 
function set_nodes_bounded_by!(
    partition::MultiLevelPartition,
    dists::Vector{Int}, 
    bounding_pts::Set{GlobalID},
    visited, 
    this_area_bounding_pts, 
    vx
)
    push!(visited, vx)
    
    if vx in bounding_pts 
        push!(this_area_bounding_pts, vx)
        return 
    end 
    
    for (ni, MLN) in partition.graph.all_neighbors_map[vx] 
        if length(MLN.finer_neighbors) == 0 && look_up_district(partition, ni) in dists
             if !(ni in visited)
                set_nodes_bounded_by!(partition, dists, bounding_pts, visited, this_area_bounding_pts, ni)
            end
        end
    end
end

function lie_about_articulation_points(
    partition::MultiLevelPartition, 
    dists, 
    rng
)::Dict{GlobalID, Set{GlobalID}}
    art_pts = articulation_points(partition, dists)
    bounded_areas = nodes_bounded_by(partition, rng; dists=dists, bounding_pts=art_pts)
    biggest_area_idx = max_index(bounded_areas, x -> length(x[1])) 
    # the main area is the one with the most pcts 

    _, bounding_art_pts = bounded_areas[biggest_area_idx]

    res = Dict{GlobalID, Set{GlobalID}}()

    for art_pt in bounding_art_pts 
        bounded_areas = nodes_bounded_by(partition, rng; dists=dists, bounding_pts=Set([art_pt]))
        @assert length(bounded_areas) > 1 
        main_idx = max_index(bounded_areas, x -> length(x[1]))
        this_bicc = Set{GlobalID}()

       
        
        for i in eachindex(bounded_areas)
            if i != main_idx 
                nodes, _ = bounded_areas[i]
                union!(this_bicc, nodes)
            end
        end

        art_pt_name = node_name(partition.graph, art_pt)
        #@show art_pt_name 
        bicc_name = [node_name(partition.graph, x) for x in this_bicc] 
        #@show bicc_name

        res[art_pt] = this_bicc
    end
    return res 
end

function check(partition) 
    for i=1:partition.num_dists 
        for j=1:partition.num_dists 
            for k in partition.adjacent_fine_neighbors[i,j]
                if look_up_district(partition, k) != j 
                    @show node_name(partition.graph, k) 
                    @show i, j 
                    @show look_up_district(partition, k)
                    return false 
                end 
            end
        end
    end
    return true 
end

function make_filter_from_cycle(cycle::Vector{Int}) 
    filter = []
    for i=1:length(cycle)-1 
        push!(filter, (cycle[i], cycle[i+1]))
    end 
    push!(filter, (cycle[length(cycle)], cycle[1]))
    return filter
end 

function sample_district_cycle(
    partition::MultiLevelPartition, 
    rng::AbstractRNG 
)::Tuple{Float64, Vector{Int}}
    g = SimpleWeightedGraph(partition.adjacency)
    edges = wilson_rst(g, rng)
    
    t = Graph(edges) 
    d1, d2 = pick_two_nonedges(partition.num_dists, edges)
    path_parents = djikstra_shortest_paths(t, d1).parents 
    cycle = Vector{Int}()
    add_cycle!(cycle, d1, d2, path_parents) 

    contracted_graph = contract_vertices(g, cycle)

    t_count_contracted = exp(log_nspanning(contracted_graph))
    t_count_original = exp(log_nspanning(g))
    cycle_length = length(cycle) 
    num_nontree_edges = ne(g) - (nv(g) - 1)

    p = t_count_original * cycle_length / (t_count_contracted * num_nontree_edges)
    return cycle, p 
end

function add_cycle!(cycle::Set{Int}, d1, d2, path_parents::Vector{Int}) 
    push!(cycle, d2) 
    if d1 != d2 
        add_cycle!(cycle, d1, path_parents[d2], path_parents)
    end
end

function pick_two_nonedges(num_dists::Int, edges)
    d1 = rand(1:num_dists, rng)
    d2 = rand(1:num_dists, rng)
    if d1 == d2 || Edge(d1, d2) in edges || Edge(d2, d1) in edges 
        return pick_two_nonedges(num_dists, edges) 
    end 
    return d1, d2
end

# requires partition.graph to be a two-level graph (cty, pct)
function construct_county_split_graph(partition::MultiLevelPartition)::SimpleWeightedGraph
    adj = zeros(Int64, partition.num_dists, partition.num_dists)
    for (node_ID, dist) in partition.node_to_district
        if dist == -1 && partition.graph.ID_to_node[node_ID].node_name !== nothing && node_ID !== partition.graph.root.global_ID
            #it's a county and it's split 
            dists = Set() 
            for child_ID in get_children_IDs(partition.graph, node_ID) 
                @assert partition.node_to_district[child_ID] !== -1 
                push!(dists, partition.node_to_district[child_ID])
            end 
            @assert length(dists) >= 2 
            for i in dists 
                for j in dists 
                    if i != j 
                        adj[i,j] += 1 
                    end 
                end
            end
        end
    end

    return SimpleWeightedGraph(adj)
end

function pick_random_triangle(partition::MultiLevelPartition, rng::AbstractRNG) 
    tris = enumerate_triangles(Graph(partition.adjacency))
    tri, p = uniformly_choose(tris, rng)
    return tri, p 
end