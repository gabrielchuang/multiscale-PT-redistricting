#=function MaxSplitNodes(f) 
    return MaxSplitNodes(f) 
end =#

function initialize_constraints()
    return Dict{Type{T} where T<:AbstractConstraint, AbstractConstraint}()
end

function add_constraint!(
    constraints::Dict,
    constraint::T
) where T <: AbstractConstraint
    constraints[T] = constraint
end

function PopulationConstraint(
    graph::BaseGraph,
    num_dists::Int,
    tolerance::Float64,
    ideal_pop::Real = 0,
)::PopulationConstraint
    if ideal_pop == 0
        ideal_pop = graph.total_pop / num_dists
    end
    # no particular reason to not use floor() instead of ceil()
    min_pop = Int(ceil((1 - tolerance) * ideal_pop))
    max_pop = Int(floor((1 + tolerance) * ideal_pop))
    return PopulationConstraint(min_pop, max_pop)
end

function PopulationConstraint(
    total_pop, num_dists, tolerance
)
    ideal_pop = total_pop / num_dists
    min_pop = Int(ceil((1 - tolerance) * ideal_pop))
    max_pop = Int(floor((1 + tolerance) * ideal_pop))
    return PopulationConstraint(min_pop, max_pop)
end

function ConstrainDiscontinuousTraversals()
    return ConstrainDiscontinuousTraversals(1)
end


function satisfies_constraint(
    constraint::MaxSplitNodes, 
    partition::MultiLevelPartition, 
    node_set=[partition.graph.root.global_ID]
)
    for node_ID in node_set 
        if constraint.max_split_nodes(partition.graph, node_ID) === nothing 
            continue 
        end 
        split_children = [child_ID for child_ID in get_children_IDs(partition.graph, node_ID) 
            if partition.node_to_district[child_ID] == -1]
        if length(split_children) > constraint.max_split_nodes(partition.graph, node_ID)
            return false 
        end
        if !satisfies_constraint(constraint, partition, split_children) 
            return false 
        end
    end
    return true 
end

function satisfies_constraint(
    constraint::MaxSplitNodesByLevel, 
    partition::MultiLevelPartition
)
    for d = 1:maximum(values(partition.graph.distance_from_root))
        nodes_at_this_level = filter(x -> last(x) == d, partition.graph.distance_from_root)
        split_nodes_at_this_level = [k for (k,v) in nodes_at_this_level 
            if haskey(partition.node_to_district, k) && partition.node_to_district[k] == -1]

        if length(split_nodes_at_this_level) >= constraint.max_split_nodes 
            return false 
        end
    end
    return true
end

function satisfies_constraint(
    constraint::ConstrainDiscontinuousTraversals,
    subgraph::MultiLevelGraph, # the merged area of districts we're recom-ing 
    node_set=[subgraph.root.global_ID]
)
    @assert constraint.max_line_segments == 1 # not built yet otherwise
    #println("checking constraints on ", node_set)
    for node_ID in node_set
        if is_intact(subgraph, node_ID)
            #println(node_ID, " is intact")
            continue
        end
        cg, vmap = child_graph(subgraph, node_ID)
        if !is_connected_bf(cg)
            # NB: I am not sure if this is working. I believe that 
            # when the child graph is not connected, it returns a subset of the connected components 
            # and so this check does not succeed. 
            #println("disc")
            #@show cg, vmap
            #@show get_children_IDs(subgraph, node_ID)
            #write_graph_nodes_to_file(subgraph, "viz/disconnected.csv", depth=2)
            return false
        end
        
        if nv(cg) != length(get_children_IDs(subgraph, node_ID))
            # NB: I'm not sure if this covers all the cases, frankly, but 
            # it seems to work for now (as a proxy for connectivity checking). 
            #println("mismatch: ", nv(cg), " ", length(get_children_IDs(subgraph, node_ID)), 
            #    " ", node_ID, " ", node_name(subgraph, node_ID))
            #println("mismatch")
            #@show cg, vmap
            #@show get_children_IDs(subgraph, node_ID)

            #write_graph_nodes_to_file(subgraph, "viz/mismatch.csv", depth=subgraph.distance_from_root[node_ID])
            return false 
        end
    
        if !satisfies_constraint(constraint, subgraph, get_children_IDs(subgraph, node_ID))
            return false
        end
    end
    return true
end


function satisfies_constraint(
    constraint::PopulationConstraint,
    partition::MultiLevelPartition,
    districts::Vector{Int}=1:partition.num_dists
)
    for di in districts
        pop = partition.dist_populations[di]
        if pop < constraint.min_pop || pop > constraint.max_pop
            return false
        end
    end
    return true
end  

function satisfies_constraints(
    partition::MultiLevelPartition,
    constraints::Dict,
    update::Union{Tuple, Nothing}=nothing
)
    snapshot = nothing
    if update !== nothing
        snapshot = take_snapshot_and_update!(partition, update...)
    end

    satisfies_constraints = true 
    changed_districts, new_dist_masks, new_dist_pops = update

    if haskey(constraints, PopulationConstraint)
        satisfies_this_constraint = true 
        if update === nothing
            satisfies_this_constraint = satisfies_constraint(constraints[PopulationConstraint], partition)
        else
            satisfies_this_constraint = satisfies_constraint(constraints[PopulationConstraint], partition, changed_districts)
        end
        if !satisfies_this_constraint 
            partition.extensions[REJECTION_COUNT]["pop constraint"] += 1
        end 
        satisfies_constraints &= satisfies_this_constraint 
    end

    if haskey(constraints, ConstrainDiscontinuousTraversals)
        satisfies_this_constraint = true 
        for i in 1:partition.num_dists 
            satisfies_this_constraint &= satisfies_constraint(constraints[ConstrainDiscontinuousTraversals], partition.subgraphs[i])
        end 
        if !satisfies_this_constraint 
            partition.extensions[REJECTION_COUNT]["continuous traversal constraint"] += 1
        end 
        satisfies_constraints &= satisfies_this_constraint 
    end

    if haskey(constraints, MaxSplitNodes)
        satisfies_this_constraint = satisfies_constraint(constraints[MaxSplitNodes], partition)
        if !satisfies_this_constraint 
            partition.extensions[REJECTION_COUNT]["max split nodes"] += 1
        end 
        satisfies_constraints &= satisfies_this_constraint 
    end

    if haskey(constraints, MaxSplitNodesByLevel)
        satisfies_this_constraint &= satisfies_constraint(constraints[MaxSplitNodesByLevel], partition)
        if !satisfies_this_constraint 
            partition.extensions[REJECTION_COUNT]["max split nodes by level"] += 1
        end 
        satisfies_constraints &= satisfies_this_constraint 
    end

    if update !== nothing
        restore_from_snapshot!(partition, snapshot...) 
    end

    return satisfies_constraints
end
