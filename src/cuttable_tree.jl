struct CuttableTree
    tree::SimpleWeightedGraph
    vmap::Vector{GlobalID}
    rev_vmap::Dict{LocalID, GlobalID}
    enter::Vector{Int}
    exit::Vector{Int}
    parent::Vector{Int}
    #name_to_index::Dict{Tuple{Vararg{String}}, Int}
        # hope this isn't needed cuz its kind of ..,.. just globalid to localid 
    pop_weight::Vector{Int}
    cuttable_edges::Vector{UndirectedEdge{GlobalID}}
    cuttable_nodes::Vector{GlobalID}
end

struct MultiScaleCuttableTree
    cuttable_trees::Dict{GlobalID,CuttableTree}
    specified_edges::Dict{GlobalID, Dict{UndirectedEdge{LocalID}, UndirectedEdge{GlobalID}}}
end

function MultiScaleCuttableTree() 
    ctrees = Dict{GlobalID,CuttableTree}()
    spec_edges = Dict{GlobalID, Dict{UndirectedEdge{LocalID}, UndirectedEdge{GlobalID}}}()
    return MultiScaleCuttableTree(ctrees, spec_edges) 
end 

function tree_induced_subgraph(graph, edges)
    new_graph = SimpleWeightedGraph() 
    new_graph.weights = zeros(nv(graph), nv(graph)) 
    for e in edges 
        new_graph.weights[src(e), dst(e)] = graph.weights[src(e), dst(e)]
        new_graph.weights[dst(e), src(e)] = graph.weights[src(e), dst(e)]
    end 
    return new_graph, 1:nv(graph)
end



""""""
function sample_spanning_tree(
    simple_graph::SimpleWeightedGraph,
    rng::AbstractRNG
)
    @assert nv(simple_graph) != 0
    if nv(simple_graph) == 1
        subtree, vmap_t2sg = induced_subgraph(simple_graph, [1])
    else
        subtree_edges = wilson_rst(simple_graph, rng)
        subtree, vmap_t2sg = tree_induced_subgraph(simple_graph, subtree_edges)
    end
    return subtree, vmap_t2sg
end

function expand_ext_decorators!(G::MultiLevelGraph, 
    expanded_ext_decorators, ext_decorators) 
    for (k,v) in ext_decorators 
        function merge(_, node_ID)
            if haskey(expanded_ext_decorators, node_ID)
                merge!(expanded_ext_decorators[node_ID], v)
            else 
                expanded_ext_decorators[node_ID] = v
            end
        end 
        ancestors_postord_foreach(G, k, merge) 
    end
end

""""""
function construct_cuttable_tree!(
    multiscale_cuttable_tree::MultiScaleCuttableTree,
    subgraph::MultiLevelGraph,  
    constraints::Dict,
    node_ID::GlobalID=subgraph.root.global_ID; 
    edge_decorators::Dict{GlobalID, Dict{GlobalID, Int}}=Dict{GlobalID, Dict{GlobalID, Int}}(),
    ext_decorators::Dict{GlobalID, Dict{GlobalID, Int}}=Dict{GlobalID, Dict{GlobalID, Int}}(),
    rng::AbstractRNG=PCG.PCGStateOneseq(UInt64),
    balance::Tuple=(1,1), 
    only_cut_top_level_edge::Bool = false, 
    expanded_ext_decorators::Dict{GlobalID, Dict{GlobalID, Int}} = Dict{GlobalID, Dict{GlobalID, Int}}(),
)
    if length(expanded_ext_decorators) == 0 && length(ext_decorators) != 0
        # the expanded_ext_decorators are the ext decorators, but each node includes the total 
        # ext decorators of all its descendants. 
        expand_ext_decorators!(subgraph, expanded_ext_decorators, ext_decorators) 
    end 

    simple_graph, vmap_sg2g = child_graph(subgraph, node_ID)
    subtree, vmap_t2sg = sample_spanning_tree(simple_graph, rng)

    #@ensure nv(simple_graph) == length(get_children_IDs(subgraph, node_ID))
    #@ensure nv(subtree) == length(get_children_IDs(subgraph, node_ID))

    vmap_t2g = [vmap_sg2g[vmap_t2sg[ii]] for ii = 1:nv(subtree)] # global index

    cuttable_tree = find_cuttable_components(subtree, vmap_t2g, subgraph,
        constraints, edge_decorators, ext_decorators, expanded_ext_decorators, balance)
   
    multiscale_cuttable_tree.cuttable_trees[node_ID] = cuttable_tree
    multiscale_cuttable_tree.specified_edges[node_ID] = Dict{UndirectedEdge{LocalID}, UndirectedEdge{GlobalID}}()

    if !only_cut_top_level_edge
        for node::GlobalID in cuttable_tree.cuttable_nodes
            sub_edge_decorators = specify_edges!(multiscale_cuttable_tree, subgraph,
                vmap_t2g, edge_decorators, ext_decorators, node, rng)

            construct_cuttable_tree!(multiscale_cuttable_tree, subgraph,
                constraints, node; edge_decorators=sub_edge_decorators, 
                ext_decorators=ext_decorators, rng=rng, balance=balance, 
                expanded_ext_decorators=expanded_ext_decorators)
        end
    end 
end

""""""
# done converting 
function find_cuttable_components(
    subtree::SimpleWeightedGraph,
    vmap::Vector{Int},
    subgraph::MultiLevelGraph,
    constraints::Dict,
    edge_decorators::Dict{GlobalID, Dict{GlobalID, Int}},
    ext_decorators::Dict{GlobalID, Dict{GlobalID, Int}},
    expanded_ext_decorators::Dict{GlobalID, Dict{GlobalID, Int}},
    balance::Tuple{Int,Int}=(1,1)
)::CuttableTree
    total_pop = get_total_population(subgraph, ext_decorators)

    topological_sort_ind = 0
    visited = [false for _ = 1:nv(subtree)]
    parent = [-1 for _ = 1:nv(subtree)]

    enter = [-1 for _ = 1:nv(subtree)]
    exit = [-1 for _ = 1:nv(subtree)]
    pop_weight = [-1 for _ = 1:nv(subtree)]

    cuttable_edges = Vector{UndirectedEdge{GlobalID}}()
    cuttable_nodes = Vector{GlobalID}()

    stack = [1] # arbitrarily root at 1

    while length(stack) > 0
        curr::LocalID = pop!(stack);

        if !visited[curr] # first visit
            visited[curr] = true
            enter[curr] = topological_sort_ind
            topological_sort_ind += 1

            push!(stack, curr) # Revisit after children
            for nbr in neighbors(subtree, curr)
                if nbr != parent[curr]
                    parent[nbr] = curr
                    push!(stack, nbr)
                end
            end
        else # second visit
            exit[curr] = topological_sort_ind
            topological_sort_ind += 1

            set_subtree_pop!(pop_weight, parent, subtree, vmap,
                subgraph, edge_decorators, expanded_ext_decorators, curr)

            node_is_cuttable = is_cuttable_node(subtree, vmap, subgraph,
                edge_decorators, ext_decorators, expanded_ext_decorators, constraints, pop_weight,
                parent, curr, balance)
            if node_is_cuttable
                push!(cuttable_nodes, vmap[curr])
            end
            edge_is_cuttable = is_cuttable_edge(curr, pop_weight,
                total_pop, constraints, balance, vmap)

            if edge_is_cuttable && curr != 1
                curr_ID = vmap[curr]
                parent_ID = vmap[parent[curr]]
                #println("cut edg len ", length(cuttable_edges))
                push!(cuttable_edges, UndirectedEdge(curr_ID, parent_ID))
            end
        end
    end

    #println("fiianl cut edg len ", length(cuttable_edges))
    rev_vmap = reverse_vmap(vmap)

    return CuttableTree(subtree, vmap, rev_vmap, enter, exit, parent,
        pop_weight, cuttable_edges, cuttable_nodes)
end

"""
    Update curr's pop to indicate pop of entire subtree rooted here.
"""
# done converting 
function set_subtree_pop!(
    pop_weight::Vector{Int},
    parent::Vector{LocalID},
    subtree::SimpleWeightedGraph,
    vmap::Vector{GlobalID},
    G::MultiLevelGraph,
    edge_decorators::Dict{GlobalID, Dict{GlobalID, Int}},
    expanded_ext_decorators::Dict{GlobalID, Dict{GlobalID, Int}},
    node_local_ID::LocalID,
)
    node_ID = vmap[node_local_ID]

    pop_w = get_node_population(G, node_ID) + # pop at this node 
     # + pop of chilldren (filter to exclude parent in rooted tree)
        (neighbors(subtree, node_local_ID)
            |> x -> filter(nbr -> nbr != parent[node_local_ID], x)
            |> x -> sum([pop_weight[y] for y in x]))

    if haskey(edge_decorators, node_ID)
        pop_w += sum(values(edge_decorators[node_ID]))
    end
    if haskey(expanded_ext_decorators, node_ID)
        pop_w += sum(values(expanded_ext_decorators[node_ID]))
    end
    pop_weight[node_local_ID] = pop_w
end

# done converting 
function get_neighbors_pop(
    subtree::SimpleWeightedGraph,
    G::MultiLevelGraph, 
    vmap::Vector{GlobalID},
    edge_decorators::Dict{GlobalID, Dict{GlobalID, Int}},
    ext_decorators::Dict{GlobalID, Dict{GlobalID, Int}},
    expanded_ext_decorators::Dict{GlobalID, Dict{GlobalID, Int}},
    pop_weight::Vector{Int},
    parent::Vector{LocalID},
    node_local_ID::LocalID,
)::Dict{GlobalID,Int}
    total_pop = get_total_population(G, ext_decorators)

    nbr_cut_pops = Dict{GlobalID,Int}()
    node_ID = vmap[node_local_ID]
    if haskey(edge_decorators, node_ID)
        merge!(nbr_cut_pops, edge_decorators[node_ID])
    end
    if haskey(expanded_ext_decorators, node_ID)
        merge!(nbr_cut_pops, expanded_ext_decorators[node_ID])
    end
    for nbr in neighbors(subtree, node_local_ID)
        nbr_ID = vmap[nbr]
        @assert !haskey(nbr_cut_pops, nbr_ID)
        if parent[node_local_ID] == nbr
            nbr_cut_pops[nbr_ID] = total_pop - pop_weight[node_local_ID]
        else
            nbr_cut_pops[nbr_ID] = pop_weight[nbr]
        end
    end
    return nbr_cut_pops
end

# done converting 
function is_cuttable_node(
    subtree::SimpleWeightedGraph,
    vmap::Vector{GlobalID},
    G::MultiLevelGraph,
    edge_decorators::Dict{GlobalID, Dict{GlobalID, Int}},
    ext_decorators::Dict{GlobalID, Dict{GlobalID, Int}},
    expanded_ext_decorators::Dict{GlobalID, Dict{GlobalID, Int}},
    constraints::Dict,
    pop_weight::Vector{Int},
    parent::Vector{LocalID},
    node_local_ID::LocalID,
    balance::Tuple{Int,Int}
)::Bool 
    node_ID = vmap[node_local_ID]
    if length(get_children_IDs(G, node_ID)) == 0
        return false
    end

    nbr_cut_pops = get_neighbors_pop(subtree, G, vmap, edge_decorators,
        ext_decorators, expanded_ext_decorators, pop_weight, parent, node_local_ID)

    nbrs = collect(keys(nbr_cut_pops))
    sub_nbrs = view(nbrs, 1:length(nbrs)-1)
    nbr_combinations = Combinatorics.powerset(sub_nbrs)

    min_pop = constraints[PopulationConstraint].min_pop
    max_pop = constraints[PopulationConstraint].max_pop

    if any([nbr_cut_pops[nbr] > max(balance[1], balance[2])*max_pop for nbr in nbrs])
        return false 
    end 

    curr_pop = get_node_population(G, node_ID)
    min_pops_1, min_pops_2 = 0, 0
    min1, min2, max1, max2 = false, false, false, false
    cuttable = false
    tot_nbr_pop = sum([nbr_cut_pops[n] for n in nbrs])

    for nbr_cmb in nbr_combinations
        min_pops_1 = sum([nbr_cut_pops[n] for n in nbr_cmb])
        min_pops_2 = tot_nbr_pop - min_pops_1

        min1 = (balance[1]*min_pop <= min_pops_1 + curr_pop)
        min2 = (balance[2]*min_pop <= min_pops_2 + curr_pop)
        max1 = (balance[1]*max_pop >  min_pops_1)
        max2 = (balance[2]*max_pop >  min_pops_2)
        cuttable = (min1 && max1 && min2 && max2)

        if balance[1] != balance[2]
            min1 = (balance[2]*min_pop <= min_pops_1 + curr_pop)
            min2 = (balance[1]*min_pop <= min_pops_2 + curr_pop)
            max1 = (balance[2]*max_pop >  min_pops_1)
            max2 = (balance[1]*max_pop >  min_pops_2)
            cuttable = cuttable || (min1 && max1 && min2 && max2)
        end

        if cuttable
            return true
        end
    end
    return false 
end

#done converting 
function is_cuttable_edge(
    node_local_ID::LocalID,
    pop_weight::Vector{Int},
    total_pop::Int,
    constraints::Dict,
    balance::Tuple{Int, Int}, 
    vmap 
)::Bool
    pop_1 = pop_weight[node_local_ID]
    pop_2 = total_pop - pop_1
    
    min_pop = constraints[PopulationConstraint].min_pop
    max_pop = constraints[PopulationConstraint].max_pop

    check_11 = balance[1]*min_pop <= pop_1
    check_11 = check_11 && pop_1 <= balance[1]*max_pop
    check_22 = balance[2]*min_pop <= pop_2
    check_22 = check_22 && pop_2 <= balance[2]*max_pop

    check_12 = balance[2]*min_pop <= pop_1
    check_12 = check_12 && pop_1 <= balance[2]*max_pop
    check_21 = balance[1]*min_pop <= pop_2
    check_21 = check_21 && pop_2 <= balance[1]*max_pop

    return (check_11 && check_22) || (check_12 && check_21)
end


function get_total_population(
    G::MultiLevelGraph, 
    ext_decorators::Dict{GlobalID, Dict{GlobalID, Int}},
)::Int 
    return get_node_population(G, G.root.global_ID) + 
        sum([sum(values(d)) for (_,d) in ext_decorators])
end 

function add_decorator(
    sub_decorators::Dict{GlobalID, Dict{GlobalID, Int}},
    fine_node::GlobalID,
    nbr_ID::GlobalID,
    pop::Int
)
    if haskey(sub_decorators, fine_node)
        @assert !haskey(sub_decorators[fine_node], nbr_ID) 
        sub_decorators[fine_node][nbr_ID] = pop
    else
        sub_decorators[fine_node] = Dict(nbr_ID=>pop)
    end
end

""""""
#= 
This function is a little concerning. We want to add a decorator, pointing from 
    fine_node to nbr_ID, but to do this we need to know the relevant populations. 
    This requires knowing who's the parent and who's the child in the 
    cuttable tree at node_ID's level - since the population is either 
    pop_weights[node_ID] or total - pop_weights[node_ID]. :V 
    the lack of invariants is very suspicious but I think this is correct? 
=# 
function add_decorator(
    sub_decorators::Dict{GlobalID, Dict{GlobalID, Int}},
    cuttable_tree::CuttableTree,
    subgraph::MultiLevelGraph,
    ext_decorators::Dict{GlobalID, Dict{GlobalID, Int}},
    fine_node::GlobalID,
    node_ID::GlobalID,
    nbr_ID::GlobalID, 
    sibling_local_ID::LocalID
)
    pop = 0
    node_local_ID = cuttable_tree.rev_vmap[node_ID]

    if sibling_local_ID == cuttable_tree.parent[node_local_ID]
        pop = get_total_population(subgraph, ext_decorators)
        pop -= cuttable_tree.pop_weight[node_local_ID]
    else
        @ensure node_local_ID == cuttable_tree.parent[sibling_local_ID]
        pop = cuttable_tree.pop_weight[sibling_local_ID]
    end
    add_decorator(sub_decorators, fine_node, nbr_ID, pop)
end



""""""
# we're zooming in on node node_ID, so we need to specify the edges from node_ID 
# to its siblings. updates specified edges and returns the edge decorators associated with node_ID.
# done  
function specify_edges!(
    multiscale_cuttable_tree::MultiScaleCuttableTree,
    subgraph::MultiLevelGraph,
    vmap_t2sg::Vector{Int},
    edge_decorators::Dict{GlobalID, Dict{GlobalID, Int}},
        # when we zoom in on a node, there are some edge decs that go to coarser lvl nodes 
    exterior_decorators::Dict{GlobalID, Dict{GlobalID, Int}},
        # when we have 2 districts, we have exterior decs for each 
    node_ID::GlobalID,
    rng::AbstractRNG
)
    #if node_ID == 2719 
    #    @show edge_decorators 
    #end 

    sub_decorators = Dict{GlobalID, Dict{GlobalID, Int}}()
    parent_ID = get_parent_ID(subgraph, node_ID)
    cuttable_tree = multiscale_cuttable_tree.cuttable_trees[parent_ID]
    specified_edges = multiscale_cuttable_tree.specified_edges[parent_ID]
        #::Dict{UndirectedEdge{LocalID}, UndirectedEdge{GlobalID}} 

    @ensure node_ID ∈ vmap_t2sg 
    node_local_ID = cuttable_tree.rev_vmap[node_ID]
    
    for sibling_local_ID::LocalID in neighbors(cuttable_tree.tree, node_local_ID)
        coarser_specified_edges = multiscale_cuttable_tree.specified_edges[get_parent_ID(subgraph, node_ID)]
        this_edge = UndirectedEdge(sibling_local_ID, node_local_ID)
        if haskey(coarser_specified_edges, this_edge)
            #if, at the level above, node_local_ID and sibling_local_ID is an already-specified edge 
            # then our neighbor's global ID is specified (i.e. finer level), and nbr_global_ID is 
            # finer level. This happens if sibling was already refined.  
            nbr_ID = coarser_specified_edges[this_edge].id1 == node_ID ?
                coarser_specified_edges[this_edge].id2 : coarser_specified_edges[this_edge].id1 
        else 
            # otherwise it's just the global ID of the sibling 
            nbr_ID = vmap_t2sg[sibling_local_ID] 
        end 

        fine_node::GlobalID = refine_node(subgraph, node_ID, nbr_ID, rng)

        specified_edges[UndirectedEdge(node_local_ID, sibling_local_ID)] = UndirectedEdge(fine_node, nbr_ID)

        add_decorator(sub_decorators, cuttable_tree, subgraph, exterior_decorators,
            fine_node, node_ID, nbr_ID, sibling_local_ID)
    end

    # the node that we're expanding has decorator(s). Must "refine" the decorators  
    # and their specified edges. 
    if haskey(edge_decorators, node_ID) 
        for (nbr_ID, pop) in edge_decorators[node_ID]
            fine_node = refine_node(subgraph, node_ID, nbr_ID, rng)
            LCA = least_common_ancestor(subgraph, node_ID, nbr_ID) 
            @ensure haskey(multiscale_cuttable_tree.specified_edges, LCA)
            #@assert hasvalue(multiscale_cuttable_tree.specified_edges[node_with_this_edge_specified], UndirectedEdge(node_ID, nbr_ID))
           
            specified_edge = multiscale_cuttable_tree.specified_edges[LCA]
            p((_,v)) = v == UndirectedEdge(node_ID, nbr_ID)
            coarse_edge, _ = select(p, collect(specified_edge))
            specified_edge[coarse_edge] = UndirectedEdge(fine_node, nbr_ID)
            add_decorator(sub_decorators, fine_node, nbr_ID, pop)
        end
    end 

    # I'm pretty sure this is not needed, because we expand the ext decorators
    #=if haskey(exterior_decorators, node_ID) 
        @show "hey"
        @assert false 
        for (nbr_ID, pop) in exterior_decorators[node_ID]
            fine_node = refine_node(subgraph, node_ID, nbr_ID, rng)
            LCA = least_common_ancestor(subgraph, node_ID, nbr_ID) 
            @ensure haskey(multiscale_cuttable_tree.specified_edges, LCA)
            #@assert hasvalue(multiscale_cuttable_tree.specified_edges[node_with_this_edge_specified], UndirectedEdge(node_ID, nbr_ID))
           
            specified_edge = multiscale_cuttable_tree.specified_edges[LCA]
            p((_,v)) = v == UndirectedEdge(node_ID, nbr_ID)
            coarse_edge, _ = select(p, collect(specified_edge))
            specified_edge[coarse_edge] = UndirectedEdge(fine_node, nbr_ID)
            add_decorator(sub_decorators, fine_node, nbr_ID, pop)
        end
    end =#

    return sub_decorators
end

""""""
# pick a random child of node_ID that is a neighbor of nbr_ID. 
# "uniformly at random", weighted by edge weight 
function refine_node(
    subgraph::MultiLevelGraph,
    node_ID::GlobalID,
    nbr_ID::GlobalID,
    rng::AbstractRNG,
)::GlobalID
    edge_weight = get_edge_weight(subgraph, node_ID, nbr_ID)
    rnd_num = rand(rng)*edge_weight
    cum_edge_weight = 0

    # NB: these are the children of node_ID that border nbr_ID
    border_children = get_children_that_border(subgraph, node_ID, nbr_ID)
    for fine_node::GlobalID in border_children
        @ensure is_in_graph(subgraph, fine_node)
        cum_edge_weight += get_edge_weight(subgraph, fine_node, nbr_ID)
        if cum_edge_weight > rnd_num
            return fine_node
        end
    end

    #@show [node_name(subgraph, x) for x in border_children] 
    #@show border_children 
    @show rnd_num 
    @show cum_edge_weight
    @show edge_weight
    @show node_ID#, node_name(subgraph, node_ID)  
    @show nbr_ID#, node_name(subgraph, nbr_ID)
    #@show [get_edge_weight(subgraph, x, nbr_ID) for x in border_children]
    #@show get_finest_level_neighbors(subgraph, nbr_ID) 
    #@show get_children_IDs(subgraph, node_ID)

   #print_tree(subgraph)

    flush(stdout)
    @assert false # should never reach here 
end

#done 
function cut_edge(
    multiscale_cuttable_tree::MultiScaleCuttableTree,
    subgraph::MultiLevelGraph,
    rng::AbstractRNG; 
    restrict_node_cuts_to::Union{Nothing, Vector{Int}}=nothing
)
    if restrict_node_cuts_to === nothing 
        cuttable_edges = all_cuttable_edges(multiscale_cuttable_tree)
    else 
        cuttable_edges = restricted_cuttable_edges(multiscale_cuttable_tree, restrict_node_cuts_to)
    end 

    edge, prob_edge = uniformly_choose(cuttable_edges, rng)
    if edge === nothing
        return [], [], edge, nothing
    end
    new_dist_masks, new_dist_pops = get_cut_node_sets_w_pop(edge, multiscale_cuttable_tree, subgraph)
    return new_dist_masks, new_dist_pops, edge, prob_edge
end

#done 
function all_cuttable_edges(
    multiscale_cuttable_tree::MultiScaleCuttableTree,
)::Set{UndirectedEdge{GlobalID}}
    cuttable_edges = Set{UndirectedEdge{GlobalID}}()
    for ctree in values(multiscale_cuttable_tree.cuttable_trees)
        union!(cuttable_edges, ctree.cuttable_edges)
    end 
    return cuttable_edges 
end

function count_cuttable_edges(
    multiscale_cuttable_tree::MultiScaleCuttableTree,
)::Int
    sum([length(ctree.cuttable_edges) 
        for ctree in values(multiscale_cuttable_tree.cuttable_trees)])
end

function restricted_cuttable_edges(
    multiscale_cuttable_tree::MultiScaleCuttableTree,
    restrict_node_cuts_to::Vector{Int} 
) 
    #@show multiscale_cuttable_tree.cuttable_trees[17]
    cuttable_edges = Set{UndirectedEdge{GlobalID}}()
    for ctree in restrict_node_cuts_to
        union!(cuttable_edges, multiscale_cuttable_tree.cuttable_trees[ctree].cuttable_edges)
    end 
    return cuttable_edges 
end

@enum Assignment RED BLUE MIXED  
# given an edge to cut and all the cuttable tree information, 
# returns the masks for each partition and their populations 
function get_cut_node_sets_w_pop(
    edge::UndirectedEdge{GlobalID},
    multiscale_cuttable_tree::MultiScaleCuttableTree,
    subgraph::MultiLevelGraph,
    ext_decorators::Dict{GlobalID, Dict{GlobalID, Int}}=Dict{GlobalID, Dict{GlobalID, Int}}(),
)::Tuple{Vector{Mask}, Vector{Int}}
    @ensure depth(subgraph, edge.id1) == depth(subgraph, edge.id2)
    @ensure get_parent_ID(subgraph, edge.id1) == get_parent_ID(subgraph, edge.id2) 

    n1_ID = edge.id1 
    n2_ID = edge.id2 
    parent_ID = get_parent_ID(subgraph, n1_ID)
    cuttable_tree = multiscale_cuttable_tree.cuttable_trees[parent_ID]
    n1_local_ID = indexin(n1_ID, cuttable_tree.vmap)
    n2_local_ID = indexin(n2_ID, cuttable_tree.vmap)

    if cuttable_tree.parent[n1_local_ID] != n2_local_ID
        (n1_ID, n2_ID) = (n2_ID, n1_ID)
        (n1_local_ID, n2_local_ID) = (n2_local_ID, n1_local_ID)
    end 
  
    assignment = Dict{GlobalID, Assignment}()
    assignment[parent_ID] = MIXED
  
    n1_enter = cuttable_tree.enter[n1_local_ID][1]
    n1_exit = cuttable_tree.exit[n1_local_ID][1]


    pop1 = cuttable_tree.pop_weight[n1_local_ID][1]
    pop2 = get_total_population(subgraph, ext_decorators) - pop1


    # assign the finest-level nodes, using enter/exit times 
    for n_local_ID = 1:nv(cuttable_tree.tree)
        n_enter = cuttable_tree.enter[n_local_ID]
        n_exit = cuttable_tree.exit[n_local_ID]
        n_ID = cuttable_tree.vmap[n_local_ID]
        if n_enter >= n1_enter && n_exit <= n1_exit
            assignment[n_ID] = RED 
        else
            assignment[n_ID] = BLUE
        end
    end


    function floodfill_assign(cuttable_tree, curr::LocalID, color)
        for x in neighbors(cuttable_tree.tree, curr) 
            if !haskey(assignment, cuttable_tree.vmap[x]) 
                assignment[cuttable_tree.vmap[x]] = color 
                floodfill_assign(cuttable_tree, x, color) 
            end 
        end 
    end 

    while parent_ID != subgraph.root.global_ID
        refined_node = parent_ID 
        parent_ID = get_parent_ID(subgraph, parent_ID)
        cuttable_tree = multiscale_cuttable_tree.cuttable_trees[parent_ID]
        vmap = cuttable_tree.vmap 
        spec_edges = multiscale_cuttable_tree.specified_edges[parent_ID]

        assignment[parent_ID] = assignment[refined_node]
        for spec_edge::UndirectedEdge{LocalID} in keys(spec_edges)
            finer_spec_edge::UndirectedEdge{GlobalID} = spec_edges[spec_edge] 

            if !(refined_node == vmap[spec_edge.id1] || refined_node == vmap[spec_edge.id2])
                continue 
            end 

            # every specified edge has one end being a node that got further refined 
            @ensure refined_node == vmap[spec_edge.id1] || refined_node == vmap[spec_edge.id2]
            # this is the *other* end of the specified edge 
            outside_node_local::LocalID = refined_node == vmap[spec_edge.id1] ? spec_edge.id2 : spec_edge.id1
            outside_node::GlobalID = vmap[outside_node_local]
            @ensure outside_node == ancestor_at_depth(subgraph, finer_spec_edge.id1, depth(subgraph, outside_node)) || 
                outside_node == ancestor_at_depth(subgraph, finer_spec_edge.id2, depth(subgraph, outside_node))
            
            # the finer end of the specified edge. Child of refined_node. may be deeper than 1 level down, though. 
            # the coarser end of the specified edge needs to be assigned to be the same color as the finer end 
            finer_node = is_descendant(subgraph, finer_spec_edge.id1, outside_node) ? finer_spec_edge.id2 : finer_spec_edge.id1
            

            @ensure haskey(assignment, refined_node)
            
            refined_node_child = closest_ancestor_with(subgraph, finer_node, 
                x -> haskey(assignment, x) && (assignment[x] == RED || assignment[x] == BLUE))

            @ensure refined_node_child != 0 

            # for example: 
            # spec_edge = county1 - county2 
            # finer_spec_edge = precint2 in count1 - cblock 3 
            # refined_node = county2
            # outside_node = county1 
            # finer_node = cblock 3 
            # refined_node_child = precinct (ancestor of cblock 3, child of refined_node)

            @ensure haskey(assignment, refined_node_child) 
            @ensure assignment[refined_node_child] == RED || assignment[refined_node_child] == BLUE 
            assignment[outside_node] = assignment[refined_node_child] 
            for nbr in neighbors(cuttable_tree.tree, outside_node_local)
                if vmap[nbr] != refined_node 
                    assignment[vmap[nbr]] = assignment[refined_node_child]
                    floodfill_assign(cuttable_tree, nbr, assignment[refined_node_child])
                end 
            end 
        end 

        for i in vmap 
            @assert haskey(assignment, i) 
        end 

    end
    
    function populate_mask(
        mask::Mask, 
        desired_color::Assignment, 
        node_ID::GlobalID
    )
        if !haskey(assignment, node_ID)  
            return false 
        elseif assignment[node_ID] == desired_color 
            mask[node_ID] = nothing
            return true  
        elseif assignment[node_ID] == MIXED 
            mask[node_ID] = Set() 
            for child_ID in get_children_IDs(subgraph, node_ID)
                res = populate_mask(mask, desired_color, child_ID) 
                if res 
                    push!(mask[node_ID], child_ID) 
                end 
            end 
            return true 
        else 
            @ensure assignment[node_ID] != desired_color
            return false 
        end 
    end 

    node_set_1 = Mask() 
    populate_mask(node_set_1, RED, subgraph.root.global_ID)
    district_1 = mask_intersection(node_set_1, subgraph.mask, subgraph.parent)

    node_set_2 = Mask() 
    populate_mask(node_set_2, BLUE, subgraph.root.global_ID)
    district_2 = mask_intersection(node_set_2, subgraph.mask, subgraph.parent)

    if pop1 <= pop2 
        return ([district_1, district_2], [pop1, pop2])
    else 
        return ([district_2, district_1], [pop2, pop1])
    end
end
