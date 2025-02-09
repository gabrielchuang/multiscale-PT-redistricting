const GlobalID = Int 
const LocalID = Int 

using Base.Iterators

struct EdgeAttributes 
    connections::Int 
    length::Float64
    weight::Float64 
    # makes it easy to add new attrs hopefully
end 
function EdgeAttributes(c, l) 
    EdgeAttributes(c, Float64(l), Float64(c)) # by default, the weight is the num connections (?)
end 

struct UndirectedEdge{ID}
    # we make Edge parametric so we can have both GlobalID and LocalID edges
    id1::ID
    id2::ID
    UndirectedEdge(id1::ID,id2::ID) where {ID} = id1 > id2 ? new{ID}(id2,id1) : new{ID}(id1,id2)
end

struct MultiLevelNeighbor
    neighbor_of::GlobalID
    global_ID::GlobalID
    #edge_weight::Real # mixed/unmixed nbr weights
    edge_attributes::EdgeAttributes
    finer_neighbors::Dict{GlobalID, MultiLevelNeighbor}
end

struct MultiLevelNode
    global_ID::GlobalID #should match the ID in the base graph, if it's a leaf 
    #parent_ID::GlobalID 
    #      Ref, so that it's mutable. 
    # unfortunately this apparently makes it super slow, so let's try making it not a ref?
    children_IDs::Vector{GlobalID} #could also have Vector{MultiLevelNode} instead, no preference atm
    node_name::Union{Tuple{String, String},Nothing} #e.g. "County", "Durham" 
end

const Mask = Dict{GlobalID, Union{Nothing, Set{GlobalID}}} 

struct MultiLevelGraph #<: AbstractGraph
    root::MultiLevelNode
    mask::Mask  
    parent::Union{Nothing, MultiLevelGraph}
    node_parents::Dict{GlobalID, GlobalID}
    ID_to_node::Dict{GlobalID, MultiLevelNode}
    base_graph::BaseGraph

    fine_neighbors::Dict{GlobalID, Vector{MultiLevelNeighbor}}
        # maps global IDs to the list of their MultiLevelNeighbor (trees). 
    all_neighbors_map::Dict{GlobalID, Dict{GlobalID, MultiLevelNeighbor}}
        # allows easy access to *all* multilevelneighbors 
        # note: is only populated as fine_neighbors are needed. 
    finest_level_neighbors::Dict{GlobalID, Dict{GlobalID, EdgeAttributes}} 
        # maps global IDs to, e.g. the census blocks that they are adjacent to. 
        # if we re-quotient nodes A and B, the children and parents of A and B will not have finest_level_neighbors change at all, and the new finest_level_neighbors for A and B ought to be easily computable from their children. 
    child_graphs::Dict{GlobalID, Tuple{SimpleWeightedGraph, Vector{GlobalID}}}
        # map a node to its child graph + vmap 
    node_attributes::Dict{GlobalID, Dict{String, Any}}
    log_tree_counts::Dict{GlobalID, Real} 

    distance_from_root::Dict{GlobalID, Int} 
    # NB: ONLY use this for telling what the **relative** levels of two nodes are, never for its absolute value. 

    extras::Tuple{Dict{GlobalID, Float64}, Dict{GlobalID, Float64}, Dict{GlobalID, Int}}
end

#============ function signatures. (for easy reference) ============
is_intact(G, node or ID)::Bool 
is_in_graph(G, node or ID)::Bool
get_children_IDs(G, node or ID)::Vector{GlobalID}
get_node_population(G, node_ID)::Real 
least_common_ancestor(G, n1_ID, n2_ID)::GlobalID

get_parent_ID(G, node_ID)::GlobalID
is_descendant(G, node_ID, parent_ID)::Bool

preord_foreach(G, node_ID, f)

ancestor_at_depth(G, node_ID, depth::Int)::GlobalID
closest_ancestor_with(G, node_ID, p)::GlobalID 

take_subgraph(G, mask)::MultiLevelGraph 

node_attributes(G, node_ID)::Dict{String, Any}
get_node_population(G, node_ID)::Real 
child_graph(G, node_ID)::Tuple{SimpleWeightedGraph, Vector{GlobalID}}
    vmap is LocalID -> GlobalID
get_finest_level_neighbors(G, node_ID)::Dict{GlobalID, EdgeAttributes}
fine_neighbors(G, node_ID)::Vector{MultiLevelNeighbor}

get_children_that_border(G, node_ID, sibling_ID)::Vector{GlobalID}

get_edge_weight(G, UndirectedEdge{GlobalID})::Real 
get_edge_weight(G, node_ID1, node_ID2)::Real 

mask_preord_foreach(mask, node_ID, f) 
mask_preord_foreach_leaf_only(mask, node_ID, f) 
mask_intersection(m1, m2, G)::Mask
mask_union(m1, m2, G)::Mask
=# 


# ================= PARENTS/CHILDREN GETTERS ================= 

function get_parent_ID(G::MultiLevelGraph, node_ID::GlobalID) 
    G.node_parents[node_ID]
end 
function get_parent_ID(G::MultiLevelGraph, node::MultiLevelNode) 
    G.node_parents[node.global_ID]
end 

function get_children_IDs(G::MultiLevelGraph, node::MultiLevelNode)::Vector{GlobalID}
    @ensure is_in_graph(G, node) 
    if haskey(G.mask, node.global_ID)
        if G.mask[node.global_ID] === nothing 
            node.children_IDs
        else collect(G.mask[node.global_ID])
        end 
    else node.children_IDs
    end 
end 

function get_children_IDs(G::MultiLevelGraph, node_ID::GlobalID)::Vector{GlobalID}
    get_children_IDs(G, G.ID_to_node[node_ID])
end 

# ================= MEMBERSHIP CHECKS =================

# WARNING: this does NOT account for lower-level children of this node being split, 
# and thus is actually kind of misleading. Potential correctness issue? 
function is_intact(G::MultiLevelGraph, node::MultiLevelNode)::Bool
    if G.parent === nothing
        return true 
    elseif haskey(G.mask, node.global_ID)
        return G.mask[node.global_ID] === nothing || 
            (haskey(G.parent.mask, node.global_ID) && 
                # the parent mask has this node, and 
            G.mask[node.global_ID] == G.parent.mask[node.global_ID] && 
                # the parent mask and this mask match for this node, and 
            mask_subset_equal(G.mask, G.parent.mask, node.global_ID)
                # all children of this node are also intact??
            )
    elseif get_parent_ID(G, node) == 0 
        return false  
    else is_intact(G, G.ID_to_node[get_parent_ID(G, node)])
    end 
end

function mask_subset_equal(mask1, mask2, node) 
    if mask1[node] == mask2[node] 
        return (mask1[node] === nothing || 
            all([mask_subset_equal(mask1, mask2, x) for x in mask1[node]]))
    else 
        return false 
    end
end

function is_intact(G::MultiLevelGraph, node_ID::GlobalID)::Bool 
    is_intact(G, G.ID_to_node[node_ID])
end 

function is_in_graph(G::MultiLevelGraph, node::MultiLevelNode)::Bool
    if haskey(G.mask, node.global_ID)
        return true 
    elseif get_parent_ID(G, node) == 0 
        return false 
    else is_intact(G, get_parent_ID(G, node))
    end 
end 

function is_in_graph(G::MultiLevelGraph, node_ID::GlobalID)::Bool 
    is_in_graph(G, G.ID_to_node[node_ID])
end 

# ================= ANCESTORS, DESCENDANTS  =================

function is_descendant(G::MultiLevelGraph, child::GlobalID, parent::GlobalID)::Bool
    if child == parent
        return true 
    elseif child == 0 
        return false
    else 
        return is_descendant(G, get_parent_ID(G, child), parent) 
    end 
end 

function ancestor_at_depth(G::MultiLevelGraph, node_ID::GlobalID, desired_depth::Int)::Union{Nothing, GlobalID}
    if depth(G, node_ID) < desired_depth
        @show "warning: no ancestor at this depth"
        return nothing 
    elseif G.distance_from_root[node_ID] == desired_depth 
        return node_ID 
    else 
        return ancestor_at_depth(G, get_parent_ID(G, node_ID), desired_depth)
    end 
end 

function closest_ancestor_with(G::MultiLevelGraph, node_ID::GlobalID, p::Function)::GlobalID 
    if p(node_ID) || node_ID == 0 
        node_ID 
    else 
        closest_ancestor_with(G, get_parent_ID(G, node_ID), p) 
    end 
end 

function is_descendant_of_intact(m::Mask, G::MultiLevelGraph, node_ID::GlobalID)
    if haskey(m, node_ID) 
        return m[node_ID] === nothing
    end
    return is_descendant_of_intact(m, G, get_parent_ID(G, node_ID))   
end 

function least_common_ancestor(G::MultiLevelGraph, n1::GlobalID, n2::GlobalID) 
    for d = min(depth(G, n1), depth(G, n2)):-1:0
        if ancestor_at_depth(G, n1, d) == ancestor_at_depth(G, n2, d)
            return ancestor_at_depth(G, n1, d)
        end 
    end 
    @assert false # should not reach here 
end 

function add_ancestors!(ancestors_so_far::Vector{GlobalID}, G::MultiLevelGraph, node_ID::GlobalID)
    push!(ancestors_so_far, node_ID)
    if node_ID == G.root.global_ID 
        return  
    end 
    add_ancestors!(ancestors_so_far, G, get_parent_ID(G, node_ID))
end 

function ancestors(G::MultiLevelGraph, node_ID::GlobalID)::Vector{GlobalID}
    ancs = [] 
    add_ancestors!(ancs, G, node_ID)
    return ancs 
end 


# ================= ITERATOR UTILITIES (looping through graphs, masks, etc.) ====================

function preord_foreach(G::MultiLevelGraph, node::MultiLevelNode, f) #pre-order traversal 
    f(G,node)
    foreach(x -> preord_foreach(G, G.ID_to_node[x], f), get_children_IDs(G, node))
end 
function preord_foreach(G::MultiLevelGraph, node_ID::GlobalID, f)
    preord_foreach(G, G.ID_to_node[node_ID], f)
end

function postord_foreach(G::MultiLevelGraph, node::MultiLevelNode, f) #post-order traversal 
    foreach(x -> postord_foreach(G, G.ID_to_node[x], f), get_children_IDs(G, node))
    f(G,node)
end 

function mask_preord_foreach(mask::Mask, node_ID::GlobalID, f::Function) 
    f(node_ID) 
    if mask[node_ID] !== nothing 
        foreach(x -> mask_preord_foreach(mask, x, f), mask[node_ID])
    end 
end 

function mask_postord_foreach(mask::Mask, node_ID::GlobalID, f::Function) 
    if mask[node_ID] !== nothing  
        foreach(x -> mask_postord_foreach(mask, x, f), mask[node_ID])
    end 
    f(node_ID) 
end 

function mask_preord_foreach_leaf_only(mask::Mask, node_ID::GlobalID, f)
    f_prime(x) = (mask[x] === nothing || length(mask[x]) == 0) ? f(x) : nothing 
    mask_preord_foreach(mask, node_ID, f_prime) 
end 

function ancestors_postord_foreach(G::MultiLevelGraph, node_ID::GlobalID, f) 
    f(G, node_ID) 
    if node_ID !== G.root.global_ID 
        ancestors_postord_foreach(G, G.node_parents[node_ID], f)
    end 
end 

# ================= MASK UTILITIES ====================

function clean_mask!(m::Mask, G::MultiLevelGraph) 
    #first, remove all "dead" leaf nodes (with children = []) and their parents, etc. 
    clean_empty(x) = 
        if m[x] !== nothing && length(m[x]) == 0 
            delete!(m, x) #remove this entry in the mask 
            delete!(m[get_parent_ID(G,x)], x) #remove x from it's parent in mask
            # note: if you're getting an exn thrown here, that means that your 
            # mask is Empty. if you need this to work, future me,
            # consider what you want the representaiton for an empty mask to be. Empty dict?  
        end 
    #postord ensures that children remove themselves from parents first 
    mask_postord_foreach(m, G.root.global_ID, clean_empty)

    clean_full(x) = 
        if m[x] !== nothing && 
            length(m[x]) == length(G.ID_to_node[x].children_IDs) && 
            # note: above we use literally .children_IDs rather than get_children_IDs 
            all(m[y] === nothing for y in m[x]) 
            mask_postord_foreach(m, x, child_of_x -> delete!(m, child_of_x))
            m[x] = nothing 
        end 
    #second, consolidate all "complete" child lists to nothing and remove extraneous stuff
    mask_postord_foreach(m, G.root.global_ID, clean_full) 
end 

#returns false if m[curr_node] should be [] 
function set_mask_intersection!(m::Mask, m1::Mask, m2::Mask, curr_node::GlobalID) 
    @ensure haskey(m1, curr_node) && haskey(m2, curr_node) 
    if m1[curr_node] === nothing  
        mask_preord_foreach(m2, curr_node, x->m[x] = m2[x])
    elseif m2[curr_node] === nothing 
        mask_preord_foreach(m1, curr_node, x->m[x] = m1[x])
    else 
        nodes_in_both = intersect(m1[curr_node], m2[curr_node]) 
        foreach(x -> set_mask_intersection!(m, m1, m2, x), nodes_in_both)
        m[curr_node] = nodes_in_both
    end
end 

function mask_intersection(m1::Mask, m2::Mask, G::MultiLevelGraph) 
    new_mask = Mask() 
    set_mask_intersection!(new_mask, m1, m2, G.root.global_ID)
    clean_mask!(new_mask, G) 
    return new_mask 
end 

function mask_union(m1::Mask, m2::Mask, G::MultiLevelGraph) 
    combine(c1, c2) = (c1 === nothing || c2 === nothing) ?
            nothing : union(c1, c2)
    new_mask = mergewith(combine, m1, m2)
    clean_mask!(new_mask, G)
    return new_mask 
end 

function mask_intersect(m1::Mask, m2::Mask, G::MultiLevelGraph)
    m = Mask() 
    for k in union(keys(m1), keys(m2)) 
        if haskey(m1, k) && !haskey(m2, k) 
            if is_descendant_of_intact(m2, G, k) 
                m[k] = m1[k] 
            end 
        elseif haskey(m2, k) && !haskey(m1, k)
            if is_descendant_of_intact(m1, G, k) 
                m[k] = m2[k] 
            end 
        elseif m1[k] === nothing 
            m[k] = m2[k] 
        elseif m2[k] === nothing 
            m[k] = m1[k] 
        else 
            #@assert typeof(m1[k]) :> Set && typeof(m2[k]) :> Set
            m[k] = intersect(m1[k], m2[k])
        end 
    end 
    clean_mask!(m, G)
    return m
end 

function mask_all_finest_nodes(m1::Mask) 
    return filter(x -> m1[x] === nothing, keys(m1))
end

# ================= PRINT UTILITIES =================

function print_tree(G::MultiLevelGraph; max_depth=nothing) 
    pr(G, node) = 
        if max_depth === nothing || G.distance_from_root[node.global_ID] <= max_depth
            if node.node_name === nothing
                println(" "^G.distance_from_root[node.global_ID], node.global_ID, " parent=", get_parent_ID(G, node)) 
            else 
                println(" "^G.distance_from_root[node.global_ID], node.global_ID, 
                " ", node.node_name[1], " ", node.node_name[2], " parent=",  get_parent_ID(G, node), 
                " pop=", get_node_population(G, node.global_ID))      
            end 
        end  
        flush(stdout)
        print("")
    preord_foreach(G, G.root, pr)
end 

function print_mask(m::Mask, G::MultiLevelGraph) 
    pr(node) = 
        println(" "^G.distance_from_root[node], node, 
        " ", G.ID_to_node[node].node_name[2], " pop=", get_node_population(G, node)) 
    mask_preord_foreach(m, G.root.global_ID, pr)
end 

# ================= OTHER STUFF =================

function node_name(G::MultiLevelGraph, node_ID::GlobalID)::Vector{String} 
    node = G.ID_to_node[node_ID]
    if node == G.root 
        return ["root"]
    else 
        parent_name = node_name(G, get_parent_ID(G, node))
        if node.node_name !== nothing 
            push!(parent_name, node.node_name[2])
        else 
            push!(parent_name, "g"*string(node_ID))
        end 
        return parent_name 
    end 
end 

function node_name_exclude_quotient_nodes(G::MultiLevelGraph, node_ID::GlobalID)::Vector{String}
    node = G.ID_to_node[node_ID]
    if node == G.root 
        return []
    else 
        parent_name = node_name_exclude_quotient_nodes(G, get_parent_ID(G, node))
        if node.node_name !== nothing 
            push!(parent_name, node.node_name[2])
        end 
        return parent_name 
    end 
end 

function fix_depths!(G::MultiLevelGraph) 
    set_depth(_, node) = if node != G.root 
        G.distance_from_root[node.global_ID] = G.distance_from_root[get_parent_ID(G, node)] + 1
    end 
    preord_foreach(G, G.root, set_depth) 
end  
 
# ============ Construct MultiLevelNode from BaseGraph ============
function build_multi_level_node(
    base_graph::BaseGraph, 
    nodes::Vector{Dict{String,Any}}, #think of this as a subset of base_graph.node_attributes 
    levels::Vector{String}, 
    current_level::Int, 
    node_parents::Dict{GlobalID, GlobalID},
    ID_to_node::Dict{GlobalID, MultiLevelNode}, 
    distance_from_root::Dict{GlobalID, Int},
    node_attributes::Dict{GlobalID, Dict{String, Any}}, 
    parent_ID::GlobalID, 
    ID_counter::Ref{GlobalID}, 
    depth::Int
)::MultiLevelNode
    current_level_name = levels[current_level] # e.g. "county"    
    node_name = nodes[1][current_level_name] # nodes[i]["county"] should be equal forall i
    children = Vector{GlobalID}() 

    if current_level == 1 #base case; we're at a leaf (e.g. census block)
        @ensure length(nodes) == 1
        node_ID = haskey(nodes[1], "id") ? nodes[1]["id"] : nodes[1]["index"]
        node_attributes[node_ID] = nodes[1]
    else 
        next_level_name = levels[current_level-1]
        node_ID = ID_counter[]; ID_counter[] += 1

        for smaller_unit in unique(map(x -> x[next_level_name], nodes))
            child_node = build_multi_level_node(
                base_graph, filter(x -> x[next_level_name] == smaller_unit, nodes), 
                levels, current_level-1, node_parents, ID_to_node, distance_from_root, 
                node_attributes, node_ID, ID_counter, depth+1)
            push!(children, child_node.global_ID)
        end 
    end

    node_parents[node_ID] = parent_ID 
    this_node = MultiLevelNode(node_ID, children, (current_level_name, node_name))
    ID_to_node[node_ID] = this_node 
    distance_from_root[node_ID] = depth  
    return this_node
end 

function MultiLevelGraph(
    base_graph::BaseGraph, 
    levels::Vector{String}
)
    nodes = deepcopy(base_graph.node_attributes)
    push!(levels, "root")
    foreach(x -> x["root"] = "root", nodes)
    foreach(i -> nodes[i]["index"] = i, 1:length(nodes))

    starting_node_ID = haskey(nodes[1], "id") ? maximum([x["id"] for x in nodes]) : length(nodes)

    ID_to_node = Dict{GlobalID, MultiLevelNode}() 
    distance_from_root = Dict{GlobalID, Int}() 
    node_parents = Dict{GlobalID, GlobalID}() 
    node_attributes = Dict{GlobalID, Dict{String, Any}}()

    root = build_multi_level_node(
        base_graph, vec(nodes), levels, 
        length(levels), node_parents, ID_to_node, distance_from_root, node_attributes, 
        0, Ref{Int}(starting_node_ID+1), 0)
    
    mask = Dict{GlobalID, Union{Nothing, Set{GlobalID}}}([root.global_ID => nothing])
    fine_neighbors = Dict{GlobalID, Vector{MultiLevelNeighbor}}()
    all_neighbors_map = Dict{GlobalID, MultiLevelNeighbor}()
    child_graphs = Dict{GlobalID, SimpleWeightedGraph}()
    tree_counts = Dict{GlobalID, Int}() 

    finest_level_neighbors = Dict{GlobalID, Vector{GlobalID}}() # need to compute this 
    extras= (Dict{GlobalID, Float64}(), Dict{GlobalID, Float64}(), Dict{GlobalID, Int}())

    G = MultiLevelGraph(root, mask, nothing, node_parents, ID_to_node, base_graph, fine_neighbors, 
        all_neighbors_map, finest_level_neighbors, child_graphs, 
        node_attributes, tree_counts, distance_from_root, extras) 

    compute_finest_level_neighbors!(G)
    return G 
end 

function finest_level_edges(G::MultiLevelGraph)::Set{UndirectedEdge} 
    all_edges = Set{UndirectedEdge}()
    for node_ID in keys(G.ID_to_node) 
        if is_in_graph(G, node_ID) && length(get_children_IDs(G, node_ID)) == 0
            fine_neighbors(G, node_ID)
            for (nbr_ID, y) in G.all_neighbors_map[node_ID] 
                if length(y.finer_neighbors) == 0 
                    push!(all_edges, UndirectedEdge(node_ID,nbr_ID))
                end
            end
        end
    end
    return all_edges
end

function finest_level_edges_adjlist(G::MultiLevelGraph)::Dict{GlobalID, Set{GlobalID}} 
    all_edges = Dict{GlobalID, Set{GlobalID}}([x => Set{GlobalID}() for x in finest_level_nodes(G)])
    for node_ID in keys(G.ID_to_node) 
        if is_in_graph(G, node_ID) && length(get_children_IDs(G, node_ID)) == 0
            fine_neighbors(G, node_ID)
            for (nbr_ID, y) in G.all_neighbors_map[node_ID] 
                if length(y.finer_neighbors) == 0 
                    push!(all_edges[node_ID], nbr_ID)
                end
            end
        end
    end
    return all_edges
end

function finest_level_nodes(G::MultiLevelGraph; is_flat=false)::Set{GlobalID}
    if is_flat && G.parent === nothing
        return setdiff(Set(keys(G.ID_to_node)), Set([G.root.global_ID]))
    elseif is_flat 
        return G.mask[G.root.global_ID]
    end

    all_nodes = Set{GlobalID}()
    for node_ID in keys(G.ID_to_node) 
        if is_in_graph(G, node_ID) && length(get_children_IDs(G, node_ID)) == 0
            push!(all_nodes, node_ID)
        end
    end
    return all_nodes
end

function is_finest_level_node(G::MultiLevelGraph, node_ID::GlobalID) 
    return length(get_children_IDs(G, node_ID)) == 0
end 

function MultiLevelGraph(
    base_graph::BaseGraph, 
    levels::Vector{String}, 
    hierarchy_path::String
)
    hierarchy = JSON.parsefile(hierarchy_path) 
    root_ID = hierarchy["metadata"]["root"]
    num_nodes = hierarchy["metadata"]["num_nodes"]

    ID_to_node = Dict{GlobalID, MultiLevelNode}() 

    for (node_ID_str, node_name) in hierarchy["nodes"]
        if length(node_name) == 0 
            @show node_ID_str
        end 
        level = levels[length(node_name)]
        node_ID = parse(Int64, node_ID_str)
        ID_to_node[node_ID] = MultiLevelNode(node_ID, [], (level, node_name[end]))
    end 

    for (node_ID_str, children) in hierarchy["children"] 
        node_ID = parse(Int64, node_ID_str)
        if haskey(ID_to_node, node_ID) 
            append!(ID_to_node[node_ID].children_IDs, children)
        else 
            ID_to_node[node_ID] = MultiLevelNode(node_ID, children, nothing)
        end
    end

    node_attributes = Dict{GlobalID, Dict{String, Any}}()
    for (node_ID, attr) in enumerate(base_graph.node_attributes)
        node_attributes[node_ID] = attr 
    end 
    distance_from_root = Dict{GlobalID, Int}(root_ID => 0, 0 => -1) # must fill 
    node_parents = Dict{GlobalID, GlobalID}(root_ID => 0) # must fill 
    mask = Dict{GlobalID, Union{Nothing, Set{GlobalID}}}([root_ID => nothing]) #filled 
    fine_neighbors = Dict{GlobalID, Vector{MultiLevelNeighbor}}() #starts empty 
    all_neighbors_map = Dict{GlobalID, MultiLevelNeighbor}() #starts empty 
    child_graphs = Dict{GlobalID, SimpleWeightedGraph}() # starts empty 
    tree_counts = Dict{GlobalID, Int}() # starts empty 
    finest_level_neighbors = Dict{GlobalID, Vector{GlobalID}}() # must fill 
    extras= (Dict{GlobalID, Float64}(), Dict{GlobalID, Float64}(), Dict{GlobalID, Int}())

    G = MultiLevelGraph(ID_to_node[root_ID], mask, nothing, node_parents, 
        ID_to_node, base_graph, fine_neighbors, 
        all_neighbors_map, finest_level_neighbors, child_graphs, 
        node_attributes, tree_counts, distance_from_root, extras) 

    function populate_node_parents(_, node) 
        for child_ID in node.children_IDs
            node_parents[child_ID] = node.global_ID 
        end 
    end 
    function populate_distance_from_root(_, node)
        distance_from_root[node.global_ID] = distance_from_root[node_parents[node.global_ID]] + 1
    end 
    preord_foreach(G, root_ID, populate_node_parents)
    preord_foreach(G, root_ID, populate_distance_from_root)
    compute_finest_level_neighbors!(G)
    return G     
end

# ============ subgraphing ============ #

function take_subgraph(
    G::MultiLevelGraph, 
    mask::Mask; 
    root=G.root
) 
    extras= (Dict{GlobalID, Float64}(), Dict{GlobalID, Float64}(), Dict{GlobalID, Int}())
    G = MultiLevelGraph(
        root === nothing ? G.root : root, 
        mask, 
        G, 
        G.node_parents,
        G.ID_to_node, 
        G.base_graph,
        Dict{GlobalID, Vector{MultiLevelNeighbor}}(), 
        Dict{GlobalID, MultiLevelNeighbor}(), 
        Dict{GlobalID, Vector{GlobalID}}(), 
        Dict{GlobalID, SimpleWeightedGraph}(), 
        Dict{GlobalID, Dict{String, Any}}(), 
        Dict{GlobalID, Int}(),
        G.distance_from_root, 
        extras
    )
    return G
end 

# ============ node_attributes ============ #

function merge_node_attrs(attrs)
    this_node_attrs = Dict{String, Any}() 
    for k in keys(attrs[1])
        if k == "min_lon" || k == "min_lat" 
            this_node_attrs[k] = minimum([c[k] for c in attrs])
        elseif k == "max_lon" || k == "max_lat" 
            this_node_attrs[k] = maximum([c[k] for c in attrs])    
        elseif typeof(attrs[1][k]) <: Real 
            #population, border etc. 
            this_node_attrs[k] = sum([c[k] for c in attrs])
        elseif length(attrs) >= 2 && 
                haskey(attrs[2], k) && attrs[1][k] == attrs[2][k]
                #this is e.g. a county name and we are a county 
                this_node_attrs[k] = attrs[1][k]
        end # otherwise, e.g. this is a prec_id and we are a county 
    end 
    return this_node_attrs
end

function node_attributes(G::MultiLevelGraph, node_ID::GlobalID)::Dict{String, Any}
    if haskey(G.node_attributes, node_ID) 
        return G.node_attributes[node_ID] 
    elseif G.parent !== nothing && is_intact(G, node_ID)
        G.node_attributes[node_ID] = node_attributes(G.parent, node_ID)
        return G.node_attributes[node_ID]
    else 
        child_node_attrs = map(y -> node_attributes(G, y), get_children_IDs(G, node_ID))

        @ensure length(child_node_attrs) >= 1 

        this_node_attrs = merge_node_attrs(child_node_attrs)  
        G.node_attributes[node_ID] = this_node_attrs
        return G.node_attributes[node_ID]
    end 
end 

function populate_extras!(G::MultiLevelGraph)
    d = (Dict{GlobalID, Float64}(), Dict{GlobalID, Float64}(), Dict{GlobalID, Int}())
    for n in finest_level_nodes(G, is_flat=false)
        @assert is_in_graph(G, n)
        #print(n, " ")
        d[1][n] = get_node_area(G, n)
        d[2][n] = get_node_perimeter(G, n)
        d[3][n] = get_node_population(G, n)
    end

    n = G.root.global_ID
    d[1][n] = get_node_area(G, n)
    d[2][n] = get_node_perimeter(G, n)
    d[3][n] = get_node_population(G, n)

    for n in keys(d[1])
        G.extras[1][n] = d[1][n]
        G.extras[2][n] = d[2][n]
        G.extras[3][n] = d[3][n]
    end
end

function remove_extras!(G::MultiLevelGraph) 
    for d in G.extras 
        for k in keys(d)
            delete!(d, k)
        end
    end
end

function get_node_population(
    G::MultiLevelGraph, 
    node_ID::GlobalID;
    merged_nodes = nothing
)::Float64 
    if length(G.extras[3]) != 0 
        return G.extras[3][node_ID]
    end
    if merged_nodes !== nothing && node_ID in keys(merged_nodes) 
        sum([node_attributes(G, x)[G.base_graph.pop_col] for x in merged_nodes[node_ID]])
    else
        node_attributes(G, node_ID)[G.base_graph.pop_col]
    end
end 

function get_node_area(G::MultiLevelGraph, node_ID::GlobalID; 
    merged_nodes = nothing
)::Real 
    @assert is_in_graph(G, node_ID)
    if length(G.extras[1]) != 0 
        return G.extras[1][node_ID]
    end
    if merged_nodes !== nothing && node_ID in keys(merged_nodes) 
        sum([node_attributes(G, x)[G.base_graph.area_col] for x in merged_nodes[node_ID]])
    else
        node_attributes(G, node_ID)[G.base_graph.area_col]
    end
end 

function get_node_bbox_area(G::MultiLevelGraph, node_ID::GlobalID) 
    R = 6371
    max_lon = node_attributes(G, node_ID)["max_lon"] 
    min_lon = node_attributes(G, node_ID)["min_lon"] 
    max_lat = node_attributes(G, node_ID)["max_lat"] 
    min_lat = node_attributes(G, node_ID)["min_lat"]
    return abs(sin(pi/180 * max_lat) - sin(pi/180 * min_lat)) * abs(max_lon - min_lon) * pi/180 * R^2
end

function get_node_sqbox_area(G::MultiLevelGraph, node_ID::GlobalID) 
    R = 6371
    max_lon = node_attributes(G, node_ID)["max_lon"] 
    min_lon = node_attributes(G, node_ID)["min_lon"] 
    max_lat = node_attributes(G, node_ID)["max_lat"] 
    min_lat = node_attributes(G, node_ID)["min_lat"]

    vert_del = abs(sin(pi/180 * max_lat) - sin(pi/180 * min_lat)) * R 
    horz_del =  abs(max_lon - min_lon) * pi/180 * R
    return vert_del > horz_del ? vert_del^2 : horz_del^2
end


function get_node_perimeter2(G::MultiLevelGraph, node_ID::GlobalID)::Real 
    if length(G.extras[2]) != 0 
        return G.extras[2][node_ID]
    end
    _ = fine_neighbors(G, node_ID)
    total = 0 
    for (_,v) in G.all_neighbors_map[node_ID]
        if length(v.finer_neighbors) == 0 # is a finest level node 
            total += v.edge_attributes.length 
        end
    end
    if "id" in keys(G.base_graph.node_attributes[1])
        for v in G.base_graph.node_attributes
            if v["id"] == node_ID 
                return total + v["border_length"]
            end
        end
    else 
        for (i,v) in enumerate(G.base_graph.node_attributes)
            if i == node_ID 
                return total + v["border_length"]
            end
        end
    end
    println(node_ID)
    return total 
end

function get_node_perimeter(G::MultiLevelGraph, node_ID::GlobalID) 
    if length(G.extras[2]) != 0 
        return G.extras[2][node_ID]
    end
    return sum([x.edge_attributes.length for x in fine_neighbors(G, node_ID)])
end

# ============ computation for child_graphs ============ #

## WARNING: only guaranteed to work when the child graph is fully connected 
# in particular, definitely fails when there is a singleton 
function child_graph(G::MultiLevelGraph, node_ID::GlobalID)
    if haskey(G.child_graphs, node_ID) 
        return G.child_graphs[node_ID] 
    elseif G.parent !== nothing && is_intact(G, node_ID)
        G.child_graphs[node_ID] = child_graph(G.parent, node_ID)
        return G.child_graphs[node_ID]
    else 
        @ensure is_in_graph(G, node_ID)

        child_IDs = get_children_IDs(G, node_ID)
        reverse_map = reverse_vmap(child_IDs)
        if length(child_IDs) == 0 
            G.child_graphs[node_ID] = (SimpleWeightedGraph(0), [])
            return G.child_graphs[node_ID] 
        end

        function mse2(child_ID) 
            nbrs = fine_neighbors(G,child_ID)
            nbrs2 = filter(x -> get_parent_ID(G, x.global_ID) == node_ID, nbrs)
            nbrs3 = map(x -> (x.global_ID, child_ID, x.edge_attributes), nbrs2)
            return nbrs3
        end

        function mse3(child_ID) 
            nbrs = [] 
            for (neigh_ID,v) in get_finest_level_neighbors(G,child_ID) #take all finest level neighbors
                if depth(G, neigh_ID) >= depth(G, child_ID)
                    potential_sibling = ancestor_at_depth(G, neigh_ID, G.distance_from_root[child_ID])
                    if get_parent_ID(G, potential_sibling) == node_ID
                        push!(nbrs, (potential_sibling, child_ID, v))
                    end
                end
            end
            return nbrs 
        end

        function mse4()
            rev_vmap = reverse_map
            srcs = Vector{Int64}([]) 
            dsts = Vector{Int64}([])
            wgts = Vector{Float64}([])
            for child_ID in child_IDs
                for (neigh_ID,v) in get_finest_level_neighbors(G,child_ID) #take all finest level neighbors
                    if depth(G, neigh_ID) >= depth(G, child_ID)
                        potential_sibling = ancestor_at_depth(G, neigh_ID, G.distance_from_root[child_ID])
                        if get_parent_ID(G, potential_sibling) == node_ID && potential_sibling < child_ID
                            push!(srcs, rev_vmap[potential_sibling])
                            push!(dsts, rev_vmap[child_ID])
                            push!(wgts, v.weight)
                        end
                    end
                end
            end
            if length(srcs) == 0 
                G.child_graphs[node_ID] = (SimpleWeightedGraph(length(child_IDs)), child_IDs)
            else 
                G.child_graphs[node_ID] = (SimpleWeightedGraph(srcs, dsts, wgts), child_IDs) 
            end  
        end

        mse4() 
        return G.child_graphs[node_ID]

        #I feel like mse2 should be faster, but 
        # make_sibling_edges is a lot faster. 
        # in either case both take a lot of the footprint 
        # of child_graph. 
        make_sibling_edges(child_ID) = (get_finest_level_neighbors(G,child_ID) #take all finest level neighbors
            |> x -> filter(kv -> depth(G, kv.first) >= depth(G, child_ID), x)
            |> x -> map(neigh_id -> (ancestor_at_depth(G, neigh_id, 
                G.distance_from_root[child_ID]), child_ID, x[neigh_id]), collect(keys(x))) 
                #get their ancestor at node_ID's depth +1 
            |> x -> filter(neighbor -> get_parent_ID(G, first(neighbor)) == node_ID, x)) 
                # should have parent = node_ID =#

        edges = (map(make_sibling_edges, child_IDs)
        |> x -> collect(flatten(x)) # flatten each child's edges into one list of edges 
        |> x -> map(((v1, v2, edge_attrs),)->(reverse_map[v1], reverse_map[v2], edge_attrs), x)
            # reverse_map is the vmap. map to vxs as a prefix of \mathbb{N}
        |> x -> filter(((v1, v2, edge_attrs),)->v1 < v2, x))
            # we have duplicates of each edge (one from each end); count each only once 
        
        if length(edges) == 0 
            G.child_graphs[node_ID] = (SimpleWeightedGraph(length(child_IDs)), child_IDs)
        else 
            G.child_graphs[node_ID] = (SimpleWeightedGraph(
                [x[1] for x in edges], [x[2] for x in edges], #endpoints of the edge 
                [x[3].weight for x in edges]), #the weight
                child_IDs) # the vmap     
        end 

        #if nv(G.child_graphs[node_ID][1]) != length(get_children_IDs(G, node_ID))
        #    @show G.distance_from_root[node_ID]
        #    write_graph_nodes_to_file(G, "viz/depthn.csv", depth=G.distance_from_root[node_ID]) 
        #    @assert false 
        #end
        return G.child_graphs[node_ID] 
    end 
end 


# ============ computation for fine/finest neighbors ============ #

function combine_edge_attrs(e1, e2)
    EdgeAttributes(e1.connections + e2.connections, e1.length + e2.length)
end 

function get_finest_level_neighbors(G::MultiLevelGraph, node_ID::GlobalID)
    if haskey(G.finest_level_neighbors, node_ID) 
        return G.finest_level_neighbors[node_ID] 
    elseif (G.parent !== nothing 
        && haskey(G.parent.finest_level_neighbors, node_ID) 
        && is_intact(G, node_ID))

        G.finest_level_neighbors[node_ID] = 
            filter(x -> is_in_graph(G, first(x)), G.parent.finest_level_neighbors[node_ID])
        return G.finest_level_neighbors[node_ID]
    elseif length(get_children_IDs(G, node_ID)) == 0 #leaf; inherit from parent 
        G.finest_level_neighbors[node_ID] = 
            filter(x -> is_in_graph(G, first(x)), get_finest_level_neighbors(G.parent, node_ID))
        return G.finest_level_neighbors[node_ID]
    else 
        compute_finest_neighbors!(G, G.ID_to_node[node_ID])
        return G.finest_level_neighbors[node_ID]
    end 
end 

function combine_merge(d1, d2) 
    merge(combine_edge_attrs, d1, d2)
end 

function iterate_merge(x) 
    d = Dict() 
    for di in x
        merge!(combine_edge_attrs, d, di)
    end
    return d
end 

function compute_finest_neighbors!(G, node) 
    if length(get_children_IDs(G, node)) > 0 #i.e. not a leaf 
    (get_children_IDs(G, node)
    |> (x -> (map(child -> get_finest_level_neighbors(G, child), x)))
    |> (x -> iterate_merge(x))
    |> (x -> (filter(neighbor -> !is_descendant(G, first(neighbor), node.global_ID), x)))
    |> (x -> (G.finest_level_neighbors[node.global_ID] = x)))
    end 
end 

function compute_finest_level_neighbors!(G::MultiLevelGraph)
    for id in keys(G.ID_to_node)
        G.finest_level_neighbors[id] = Dict{GlobalID, EdgeAttributes}()
    end 

    # initialize the finest level's finest_level_neighbors (which is just normal neighbors)
    for (edge,attrs) in G.base_graph.edge_attributes
        cedge = collect(edge) # potentially wastage of memory/allocation 
        G.finest_level_neighbors[cedge[2]][cedge[1]] = EdgeAttributes(attrs["connections"], attrs["length"])
        G.finest_level_neighbors[cedge[1]][cedge[2]] = EdgeAttributes(attrs["connections"], attrs["length"])
    end 
    postord_foreach(G, G.root, compute_finest_neighbors!)
end 

function construct_multi_level_neighbor(
    G::MultiLevelGraph, 
    children_map::Dict{GlobalID, Set{GlobalID}}, 
    neighbors::Dict{GlobalID, EdgeAttributes}, 
    neighbor_of::GlobalID, base::GlobalID
)
    if length(children_map[base]) == 0 
        G.all_neighbors_map[neighbor_of][base] = MultiLevelNeighbor(neighbor_of, base, neighbors[base], Dict{GlobalID, MultiLevelNeighbor}()) 
    else 
        finer_neighbors = Vector{MultiLevelNeighbor}()
        for c in children_map[base]
            nbr = construct_multi_level_neighbor(G, children_map, neighbors, neighbor_of, c)
            push!(finer_neighbors, nbr)
        end 
        finer_neighbors_dict = Dict([x.global_ID => x for x in finer_neighbors])
        edge_attr = (map(x -> x.edge_attributes, finer_neighbors) 
            |> x -> reduce(combine_edge_attrs, x))

        G.all_neighbors_map[neighbor_of][base] = MultiLevelNeighbor(neighbor_of, base, edge_attr, finer_neighbors_dict) 
    end 
    return G.all_neighbors_map[neighbor_of][base]
end 

function compute_all_fine_neighbors(G::MultiLevelGraph) 
    preord_foreach(G, G.root, (_, x) -> fine_neighbors(G, x.global_ID))
end 

# i BELIEVE that the roots of the multilevelneighbors of some node X are the 
# children of the LCA of adjacent districts. I.e., siblings, or uncle nodes 
# in the hierarchical tree. 
function fine_neighbors(G::MultiLevelGraph, node_ID::GlobalID)
    if haskey(G.fine_neighbors, node_ID) 
        @ensure haskey(G.all_neighbors_map, node_ID)
        return G.fine_neighbors[node_ID] 
    elseif G.parent !== nothing && 
        haskey(G.parent.fine_neighbors, node_ID) && 
        is_intact(G, node_ID) && 
        all([is_intact(G, x.global_ID) for x in G.parent.fine_neighbors[node_ID]]) 
            G.fine_neighbors[node_ID] = G.parent.fine_neighbors[node_ID]
            G.all_neighbors_map[node_ID] = G.parent.all_neighbors_map[node_ID]

            return G.fine_neighbors[node_ID]  

        # question: do we want to try to save computation if this node is 
        # NOT intact? Unsure. I suspect it is not worth it, since 
        # the computation of fine_neighbors is actually not too bad. 
        # (asymptotically about nd log n, where n = # finest level neighbors 
        # and d = depth of tree rooted at node_ID)
    else 
        #compute from scratch 
        finest_lvl_neighbors = filter(x -> is_in_graph(G, first(x)),
            get_finest_level_neighbors(G, node_ID))
        n_iter = collect(keys(finest_lvl_neighbors))

        neighbor_children = Dict{GlobalID, Set{GlobalID}}([x => Set{GlobalID}() for x in n_iter])
        horizontal_neighbors = Set{GlobalID}()
        for neighbor_ID in n_iter
            parent_ID = get_parent_ID(G, neighbor_ID)
            if parent_ID == 0 
                continue  
            end 
            if haskey(neighbor_children, parent_ID) 
                push!(neighbor_children[parent_ID], neighbor_ID) 
            else
                neighbor_children[parent_ID] = Set{GlobalID}(neighbor_ID)
            end
            
            if is_descendant(G, node_ID, parent_ID)
                push!(horizontal_neighbors, neighbor_ID)
                continue 
            end 
            push!(n_iter, parent_ID)
        end 

        if !haskey(G.all_neighbors_map, node_ID) 
            G.all_neighbors_map[node_ID] = Dict{GlobalID, MultiLevelNeighbor}() 
        end 
        nbrs = map(x -> construct_multi_level_neighbor(G, 
            neighbor_children, finest_lvl_neighbors, node_ID, x), collect(horizontal_neighbors)) 
        G.fine_neighbors[node_ID] = nbrs

        #=if toggle_ENSURE()
            parent_graph, vmap = child_graph(G, get_parent_ID(G, node_ID))
            node_local = indexin(node_ID, vmap)[1]
            if length(G.fine_neighbors[node_ID]) != length(neighbors(parent_graph, node_local))
                @show [node_name(G, vmap[x]) for x in neighbors(parent_graph, node_local)]
                @show [node_name(G, x.global_id) for x in G.fine_neighbors[node_ID]]
                flush(stdout)
                @assert false
            end
            @ensure length(G.fine_neighbors[node_ID]) == length(neighbors(parent_graph, node_local))
        end =#

        return G.fine_neighbors[node_ID]
    end 
end 

function depth(G::MultiLevelGraph, ID::GlobalID) 
    G.distance_from_root[ID]
end 

function get_edge_attributes(G::MultiLevelGraph, edge::UndirectedEdge{GlobalID})
    ID_high, ID_low = depth(G, edge.id1) < depth(G, edge.id2) ? 
        (edge.id1, edge.id2) : (edge.id2, edge.id1)
    _ = fine_neighbors(G, ID_high)

    # ensure that the all_neighbors_map is computed for ID_high 
    if !haskey(G.all_neighbors_map[ID_high], ID_low)
        return 0
    end 
    return G.all_neighbors_map[ID_high][ID_low].edge_attributes 
end 

function get_edge_weight(G::MultiLevelGraph, edge::UndirectedEdge{GlobalID})
    return get_edge_attributes(G, edge).weight
end 
function get_edge_connections(G::MultiLevelGraph, edge::UndirectedEdge{GlobalID})
    return get_edge_attributes(G, edge).connections
end 

function get_edge_connections(G::MultiLevelGraph, n1::GlobalID, n2::GlobalID)
    return get_edge_connections(G, UndirectedEdge(n1, n2)) 
end 
function get_edge_weight(G::MultiLevelGraph, n1::GlobalID, n2::GlobalID)
    return get_edge_weight(G, UndirectedEdge(n1, n2)) 
end 
function get_edge_attributes(G::MultiLevelGraph, n1::GlobalID, n2::GlobalID)
    return get_edge_attributes(G, UndirectedEdge(n1, n2)) 
end 

# gets child nodes that share a border with sibling_ID
function get_children_that_border(
    G::MultiLevelGraph, 
    node_ID::GlobalID, 
    sibling_ID::GlobalID
)
    if length(get_children_IDs(G, node_ID)) == 0
        return []
    else 
        return (collect(keys(get_finest_level_neighbors(G, sibling_ID)))
        |> x -> filter(nbr -> depth(G, nbr) >= depth(G, node_ID)+1, x)
        |> x -> map(nbr -> ancestor_at_depth(G, nbr, depth(G, node_ID)+1), x)
        |> x -> filter(nbr -> get_parent_ID(G, nbr) == node_ID, x))
    end 
end 

# ================== TREE COUNTS ==================
function log_tree_count(G::MultiLevelGraph, node::MultiLevelNode)::Real 
    if !haskey(G.log_tree_counts, node.global_ID) 
        if false #G.parent !== nothing && is_intact(G, node)
            @show node, "is in this case??"
            G.log_tree_counts[node.global_ID] = log_tree_count(G.parent, node)
        else
            s_graph = child_graph(G, node.global_ID)[1]
            G.log_tree_counts[node.global_ID] = 
                nv(s_graph) == 1 || !is_connected_bf(s_graph) ? 0 : log_nspanning(s_graph)
        end 
    end 
    return G.log_tree_counts[node.global_ID]
end 

function log_tree_count(G::MultiLevelGraph, node::GlobalID)::Real 
    log_tree_count(G, G.ID_to_node[node])
end 

function precompute_log_tree_counts!(G::MultiLevelGraph) 
    preord_foreach(G, G.root, log_tree_count)
end 

function log_hierarchical_tree_count(G::MultiLevelGraph)::Real 
    # requires all log tree counts to be precomputed 
    sum(values(G.log_tree_counts))
end 

function write_graph_to_file(G::MultiLevelGraph, output_file::String)
    out = []
    function add_to_out(_, node)
        node_name = node_name_exclude_quotient_nodes(G, node.global_ID)
        if length(node_name) == 2
            s = join(node_name, ",")
            push!(out, "\""*s*"\",\""*string(G.distance_from_root[node.global_ID])*"\"")
        end 
    end 
    preord_foreach(G, G.root, add_to_out)

    s = "node_name,depth\n"*join(out, "\n")
    
    open(output_file, "w") do file 
        write(file, s)
    end 
end 

function write_graph_nodes_to_file(
    G::MultiLevelGraph, 
    output_file::String; depth=3, use_global_IDs=false)
    out = []
    num_named_levels = length(node_name_exclude_quotient_nodes(
        G, G.ID_to_node[1].global_ID))

    function add_to_out(_, node)
        node_name = node_name_exclude_quotient_nodes(G, node.global_ID)
        if length(node_name) == num_named_levels && G.distance_from_root[node.global_ID] >= depth
            s = join(node_name, ",")
            output_data = ancestor_at_depth(G, node.global_ID, depth)
            push!(out, (s, output_data))
        end 
    end 
    preord_foreach(G, G.root, add_to_out)

    out_string = []
    if use_global_IDs
        for (s, data) in out 
            push!(out_string, "\""*string(s)*"\","*string(data))
        end 
    else
        remap = []
        for (s, data) in out 
            if !(data in remap)
                push!(remap, data) 
            end
            output_remap = indexin(data, remap)[1]
            push!(out_string, "\""*s*"\","*string(output_remap))
        end 
    end
    s = "node_name,node\n"*join(out_string, "\n")

    open(output_file, "w") do file 
        write(file, s)
    end 
end 

function write_hierarchy_to_file(G::MultiLevelGraph, output_file::String, metadata=Dict())
    data = Dict{String, Any}()

    dir = dirname(output_file)
    if !isdir(dir)
        mkpath(dir)
    end

    metadata["root"] = G.root.global_ID 
    metadata["num_nodes"] = length(G.ID_to_node) 

    data["metadata"] = metadata
    data["nodes"] = Dict{GlobalID, Any}() 
    data["children"] = Dict{GlobalID, Vector{GlobalID}}() 
    for i=1:length(G.ID_to_node)
        if G.ID_to_node[i].node_name !== nothing && G.ID_to_node[i] != G.root 
            data["nodes"][i] = node_name_exclude_quotient_nodes(G, i)
        end 
        data["children"][i] = G.ID_to_node[i].children_IDs
    end
    open(output_file,"w") do f
        JSON.print(f, data, 1)
    end
end

function graph_subset_by_cty(G, G_flat, incl_counties, NC_oriented_nbrs) 
    mask = Dict{GlobalID, Union{Nothing, Set{GlobalID}}}([
        G.root.global_ID => Set()])
    for node in G.root.children_IDs 
        if any([G.ID_to_node[node].node_name[2] == x for x in incl_counties])
            for child in G.ID_to_node[node].children_IDs
                push!(mask[G.root.global_ID], child)
                mask[child] = nothing 
            end
        end
    end

    Gnc3 = take_subgraph(G_flat, mask)
    @show length(finest_level_nodes(Gnc3))

    oriented_nbrs2 = Dict{GlobalID, Vector{GlobalID}}()
    for (k,v) in NC_oriented_nbrs 
        if is_in_graph(Gnc3, k) 
            oriented_nbrs2[k] = [is_in_graph(Gnc3, x) ? x : -1 for x in v]
            remove_adj_duplicates_cyclic!(oriented_nbrs2[k])
        end
    end
    return Gnc3, oriented_nbrs2
end