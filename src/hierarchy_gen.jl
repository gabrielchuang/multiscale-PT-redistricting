function generate_level(G, num_to_merge, og_onbrs, p1s, p2s, cs, pop_weight1, iso_weight)
    n = length(finest_level_nodes(G))
    pop_target = get_node_population(G, G.root.global_ID)/(n - num_to_merge) 

    #pop_weight1 = 1
    pop_weight2 = 0.01
    #iso_weight = 1

    edge_weights = Dict()
    for edge in finest_level_edges(G) 
        n1, n2 = edge.id1, edge.id2
        weight = 0

        n1pop = get_node_population(G, n1)
        n2pop = get_node_population(G, n2) 
        npop = n1pop + n2pop
        
        # a larger weight overall means higher priority 

        # this targets getting the merged node's population close to the target. we want to prioritize pairs with small values of dist_from_target. 
        dist_from_target = npop < pop_target ? pop_target/npop : npop/pop_target #always > 1; smaller means higher priority 
        weight += -pop_weight1 * dist_from_target
        push!(p1s, -pop_weight1 * dist_from_target)

        # this targets merging small nodes. We want to prioritize nodes with large values of mpop_dist_from_target. 
        mpop_dist_from_target = max(pop_target/n1pop, pop_target/n2pop) 
        weight += pop_weight2 * mpop_dist_from_target
        push!(p2s, pop_weight2 * mpop_dist_from_target)

        o_iso = 0.5 * (get_node_perimeter(G, n1)^2 / get_node_area(G, n1) + 
            get_node_perimeter(G, n2)^2 / get_node_area(G, n2))
        n_iso = merged_iso(G, n1, n2)

        if abs(o_iso - n_iso) > 10
            # this soft-caps the delta iso score. we want to penalize nodes with abs(o_iso - n_iso) > 10. 
            weight += -10 # -iso_weight * 10
        end
        iso_ratio = o_iso > n_iso ? o_iso/n_iso : n_iso/o_iso 
        # this targets merging nodes whose iso ratio remains close. We want to prioritize pairs with small values of iso_ratio. 
        weight += -iso_weight * iso_ratio
        push!(cs, -iso_weight * iso_ratio)
        
        edge_weights[edge] = weight 
    end

    sorted_by_weight = sort(collect(keys(edge_weights)), by=x -> edge_weights[x], rev=true)
    included_nodes = Set() 
    cur_onbrs = og_onbrs
    all_nodes = finest_level_nodes(G)

    to_return = [] 
    for edge in sorted_by_weight 
        if edge.id1 in included_nodes || edge.id2 in included_nodes 
            continue 
        else
            #check it's not an articulation point 
            #c_onbrs, c_pair_parent = onbrs_if_matching_made(G, [edge], og_onbrs) 
            c_onbrs, c_pair_parent = onbrs_if_matching_made_incremental(all_nodes, [edge], cur_onbrs) 
            
            # Warning: this is not exhaustive in cases with lots of corner touches. :/ 
            if count(x -> x == edge.id1, og_onbrs[edge.id2]) > 1 || 
                    count(x -> x == edge.id2, og_onbrs[edge.id1]) > 1 
                    @assert count(x -> x==edge.id1, og_onbrs[edge.id2]) > 1
                    @assert count(x -> x==edge.id2, og_onbrs[edge.id1]) > 1
                    #println("encloses a space")
                    continue
                elseif count(x -> x == -1, c_onbrs[c_pair_parent[edge.id1]]) > 1 
                    #println("touches edge twice") 
                    continue  
                else
                    push!(to_return, edge)
                    push!(included_nodes, edge.id1) 
                    push!(included_nodes, edge.id2) 
                    cur_onbrs = c_onbrs 
                    push!(all_nodes, c_pair_parent[edge.id1])
                    delete!(all_nodes, edge.id1)
                    delete!(all_nodes, edge.id2)
                end
           end

        if length(to_return) == num_to_merge
            return to_return 
        end
    end
    println("warning: returning less than requested: ", length(to_return))
    return to_return 

end

function merged_iso(G, n1, n2) 
    l1 = sum([z.length for (x,z) in get_finest_level_neighbors(G, n1)])
    l2 = sum([z.length for (x,z) in get_finest_level_neighbors(G, n2)])
    l12 = get_edge_attributes(G, UndirectedEdge(n1,n2)).length

    a1 = get_node_area(G, n1)
    a2 = get_node_area(G, n2)

    return (l1 + l2 - l12 - l12)^2 / (a1 + a2)
end

function onbrs_if_matching_made_incremental(all_nodes, matching, onbrs) 
    pair_parent = Dict{GlobalID, GlobalID}()
    new_onbrs = Dict{GlobalID, Vector{Vector{GlobalID}}}()

    matched_nodes = Set()
    for edge in matching 
        push!(matched_nodes, edge.id1) 
        push!(matched_nodes, edge.id2)
    end

    for n1 in all_nodes 
        if n1 in matched_nodes continue end 
        pair_parent[n1] = n1
        new_onbrs[n1] = [onbrs[n1]]
    end

    node_ID = maximum(all_nodes) + 1
    for edge in matching
        node_ID += 1 
        n1, n2 = edge.id1, edge.id2
        pair_parent[n1] = node_ID
        pair_parent[n2] = node_ID
        new_onbrs[node_ID] = combine_onbrs(n1, n2, onbrs)
    end

    coarsened_onbrs = Dict{GlobalID, Vector{GlobalID}}()
    for (k,v) in new_onbrs
        coarse_onbrs = []
        for section in v
            this_section_parents = [x == -1 ? -1 : pair_parent[x] for x in section]
            remove_adj_duplicates_cyclic!(this_section_parents)
            push!(coarse_onbrs, this_section_parents)
        end
        coarsened_onbrs[k] = collect(flatten(coarse_onbrs))
    end
    return coarsened_onbrs, pair_parent 
end