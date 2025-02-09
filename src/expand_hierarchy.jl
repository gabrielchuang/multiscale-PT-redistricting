function expand_hierarchy!(
    G::MultiLevelGraph,
    desired_fan_out::Int;   
    max_attempts::Int = 5, 
    rng::AbstractRNG=PCG.PCGStateOneseq(UInt64),
)
    expand(_, node) = expand_hierarchy_at_node!(G, desired_fan_out, node, max_attempts=max_attempts, rng=rng)
    preord_foreach(G, G.root, expand)
    fix_depths!(G) 
end 

function expand_hierarchy!(
    G::MultiLevelGraph,
    fan_out_predicate::Function;   
    max_attempts::Int = 5, 
    rng::AbstractRNG=PCG.PCGStateOneseq(UInt64)
)
    n = 0 
    function expand(_, node)
        n += 1 
        if mod(n, 400) == 0
            @show n
        end 
        fan_out, force_fanout = fan_out_predicate(G, node)
        expand_hierarchy_at_node!(G, fan_out, node, max_attempts=max_attempts, rng=rng, force_fanout=force_fanout)
    end 
    preord_foreach(G, G.root, expand)
    fix_depths!(G) 
end 

function expand_hierarchy_at_node!(
    G::MultiLevelGraph,
    desired_fan_out::Union{Int, Nothing}, 
        # this is a minimum on the fan-out, probably? 
    curr_node::MultiLevelNode; 
    max_attempts::Int = 5, 
    rng::AbstractRNG=PCG.PCGStateOneseq(UInt64),
    force_fanout::Bool = false, 
    pop_deviation_allowance_expansion_rate::Real = 0.3 
)
    if desired_fan_out === nothing 
        return 
    end 

    curr_node_ID = curr_node.global_ID 
    children_IDs = get_children_IDs(G, curr_node_ID)
   
    if force_fanout
        if length(children_IDs) < desired_fan_out 
            println("warning: cannot force this fanout: ", desired_fan_out)
            return 
        end 
    elseif length(children_IDs) < desired_fan_out^2
        return 
        # too few children to split into sub-nodes 
    end 

    this_node_pop = get_node_population(G, curr_node_ID)
    ideal_subnode_pop = this_node_pop / desired_fan_out
    
    # lots of q's about the optimal pop deviation. If 
    # we make it really strict, it likely won't mix well? 

    min_pop_dev = 1.05
    pop_deviation_allowance_expansion_rate = 0.1 

    for pop_deviation_allowance in min_pop_dev:pop_deviation_allowance_expansion_rate:(desired_fan_out+1)
        constraints = initialize_constraints()
        min_pop = ideal_subnode_pop * 1.0/pop_deviation_allowance
        max_pop = ideal_subnode_pop * pop_deviation_allowance
        add_constraint!(constraints, PopulationConstraint(min_pop, max_pop))

        this_node_mask = Mask() 
        this_node_mask[curr_node.global_ID] = nothing  
        this_node_subgraph = take_subgraph(G, this_node_mask; root=curr_node)

        partition = MultiLevelPartition(this_node_subgraph, constraints, desired_fan_out; 
            rng=rng, max_attempts=max_attempts, only_cut_top_level_edge=true, first_comp_take=pop_deviation_allowance < 1.3)
        if partition === nothing 
            #println("pop dev loosening from ", pop_deviation_allowance)
            #the population dev is too tight; loosen and try again 
            continue 
        end 
        new_curr_children = Set()
        for si=1:desired_fan_out
            this_subnode_subgraph::Mask = partition.district_to_nodes[si]

            subnode_children::Vector{GlobalID} = collect(this_subnode_subgraph[curr_node_ID])

            if length(subnode_children) == 1 
                # singleton sub-node; skip 
                push!(new_curr_children, subnode_children[1])
                continue 
            end 
            @assert this_subnode_subgraph[curr_node_ID] !== nothing 


            subnode_ID = maximum(keys(G.ID_to_node)) + 1 
            subnode = MultiLevelNode(subnode_ID, subnode_children, nothing)
            G.node_parents[subnode.global_ID] = curr_node_ID

            for subnode_child_ID in subnode_children 
                G.node_parents[subnode_child_ID] = subnode_ID
                #G.ID_to_node[subnode_child_ID].parent_ID[] = subnode_ID
            end 
            G.ID_to_node[subnode_ID] = subnode
            G.distance_from_root[subnode_ID] = G.distance_from_root[curr_node_ID] + 1
            fix_depths!(G)
            compute_finest_neighbors!(G, subnode)

            push!(new_curr_children, subnode_ID)
        end 

        replace!(G.ID_to_node[curr_node_ID].children_IDs, new_curr_children)
        delete!(G.child_graphs, curr_node_ID)
    
        empty!(G.fine_neighbors) 
        # this is a VERY coarse call to empty!. We should only need to delete a few 
        # entries. But expand_hierarchy should only happen once per mcmc so hopefully it is Ok 

        # (old) NOTE: we might need to clear out G.fine_neighbors, but I'm assuming that this 
        # function only ever gets called on a fresh graph. But we'll need that functionality 
        # when merge-splitting/requotienting anyway so once it's done we should add it here 
        # just in case? 
        return 
    end 
    println("warning: failed to expand partition at ", node_name(G, curr_node.global_ID)[end])
end