# note: using single node flip as it is now will kill the shared_nodes field of 
# multi level partition, so for now be careful mixing forest recom with SNF. 

# given a partially-computed energy for the pre-node-flip state, z
# compute the energy after the node flip. 
function finish_energy_with_SNF(
    initial_energies::Vector, 
    partition::MultiLevelPartition, 
    merged_nodes::Union{Nothing, Dict{GlobalID, Set{GlobalID}}},
    measure::Measure, 
    dist_to::Int, 
    pct::GlobalID
)

    dist_from = look_up_district(partition, pct) 
    if dist_from == -1 
        @assert false 
    end 
    total_energy = 0

    for i = 1:length(measure.weights)
        partial = initial_energies[i]
        if measure.scores[i].energy_type == POPULATION 
            pct_pop = get_node_population(partition.graph, pct; merged_nodes=merged_nodes)
            
            partial[dist_from] -= pct_pop 
            partial[dist_to] += pct_pop 
            total_energy += measure.weights[i] * measure.scores[i].finish_compute(partial) 
            ediff = measure.weights[i] * measure.scores[i].finish_compute(partial) 
            
            partial[dist_from] += pct_pop 
            partial[dist_to] -= pct_pop 

        elseif measure.scores[i].energy_type == COMPACTNESS
            pct_perim_dist_from = node_perim_with_district(partition, merged_nodes, pct, dist_from)
            pct_perim_dist_to = node_perim_with_district(partition, merged_nodes, pct, dist_to)
            pct_perim_total = node_perim(partition, pct; merged_nodes=merged_nodes)
            pct_area = get_node_area(partition.graph, pct; merged_nodes=merged_nodes) 

            areas, perims = partial 
            perims[dist_from] -= (pct_perim_total - 2*pct_perim_dist_from)
            perims[dist_to] += (pct_perim_total - 2*pct_perim_dist_to)
            areas[dist_from] -= pct_area 
            areas[dist_to] += pct_area 
            
            this_energy = measure.scores[i].finish_compute((areas, perims))
            ediff = measure.weights[i] * this_energy
            total_energy += ediff 

            perims[dist_from] += (pct_perim_total - 2*pct_perim_dist_from)
            perims[dist_to] -= (pct_perim_total - 2*pct_perim_dist_to)
            areas[dist_from] += pct_area 
            areas[dist_to] -= pct_area 
        else
            @assert false "no other energies implemented"
        end 
    end
    return total_energy 
end

# a speedy version when you have exactly one energy and it's compactness 
function finish_energy_with_SNF(
    partition::MultiLevelPartition, 
    weight::Float64, 
    dist_to::Int, 
    pct::GlobalID,
    og_iso::Float64
)
    dist_from = look_up_district(partition, pct) 

    areas, perims = partition.dist_areas, partition.dist_perimeters
    p_from, p_to, a_from, a_to = perims[dist_from], perims[dist_to], areas[dist_from], areas[dist_to]

    pct_perim_dist_from = node_perim_with_district(partition, nothing, pct, dist_from)
    pct_perim_dist_to = node_perim_with_district(partition, nothing, pct, dist_to)
    pct_perim_total = node_perim(partition, pct)
    pct_area = get_node_area(partition.graph, pct) 

    norm = haskey(partition.extensions, "norm") ? partition.extensions["norm"] : 1

    this_energy = og_iso^norm 
    this_energy -= (p_from^2/a_from)^norm
    this_energy -= (p_to^2/a_to)^norm
    this_energy += ((p_from - (pct_perim_total - 2*pct_perim_dist_from))^2/(a_from-pct_area))^norm
    this_energy += ((p_to + (pct_perim_total - 2*pct_perim_dist_to))^2/(a_to+pct_area))^norm 
        
    total_energy = weight * this_energy^(1/norm)

    return total_energy 
end

function articulation_points_rec(
    partition::MultiLevelPartition, 
    district::Int, 
    merged_nodes=nothing; 
    seed = nothing
)::Set{GlobalID}
    acc = Set{GlobalID}()
    visited = Set{GlobalID}() 
    depth = Dict{GlobalID, Int}()
    low = Dict{GlobalID, Int}()
    parent = Dict{GlobalID, GlobalID}()
    vx = seed === nothing ? any_pct_in(partition.subgraphs[districts[1]], partition.graph.root.global_ID) : seed  
    articulation_points_pd!(partition, districts, acc, visited, depth, low, parent, vx, 0, merged_nodes)
    return acc 
end

function any_pct_in(graph::MultiLevelGraph, node::GlobalID) 
    @ensure is_in_graph(graph, node)
    if length(get_children_IDs(graph, node)) == 0
        return node 
    else
        return any_pct_in(graph, get_children_IDs(graph, node)[1])
    end
end

# writes the articulation points to `acc`, using the biconnected componentS DFS algorithm
function articulation_points_pd!(
    partition::MultiLevelPartition, dist::Int, acc::Set{GlobalID},
    visited::Set{GlobalID}, depth::Dict{GlobalID, Int}, 
    low::Dict{GlobalID, Int}, parent::Dict{GlobalID, GlobalID},
    vx::GlobalID, d::Int,
    merged_nodes=nothing
)
    push!(visited, vx) 
    depth[vx] = d 
    low[vx] = d 
    num_children = 0 
    is_articulation = false

    for (ni, MLN) in partition.graph.all_neighbors_map[vx] 
        if length(MLN.finer_neighbors) == 0 && look_up_district(partition, ni) == dist
            if merged_nodes !== nothing && vx in keys(merged_nodes) && ni in merged_nodes[vx] 
                continue #skip the merged nodes adj to vx if vx is an art pt
            end

            if !(ni in visited)
                parent[ni] = vx 
                articulation_points_pd!(partition, dists, acc, visited, depth, low, parent, ni, d+1, merged_nodes)
                num_children += 1
                if low[ni] >= depth[vx] 
                    is_articulation = true 
                end
                low[vx] = low[vx] < low[ni] ? low[vx] : low[ni]
            elseif !haskey(parent, vx) || ni != parent[vx] 
                low[vx] = low[vx] < depth[ni] ? low[vx] : depth[ni]
            end
        end
    end
    if (haskey(parent, vx) && is_articulation) || (!haskey(parent, vx) && num_children>1)
        push!(acc, vx) 
    end
end


function articulation_points_ne(
    nodes::Set{GlobalID}, edges::Dict{GlobalID, Set{GlobalID}}
)::Set{GlobalID}
    acc = Set{GlobalID}()
    visited = Set{GlobalID}() 
    depth = Dict{GlobalID, Int}()
    low = Dict{GlobalID, Int}()
    parent = Dict{GlobalID, GlobalID}()
    vx = first(nodes)
   
    articulation_points_ne!(nodes, edges, acc, visited, depth, low, parent, vx, 0)
    return acc 
end

function articulation_points_ne!(
    nodes::Set{GlobalID}, edges::Dict{GlobalID, Set{GlobalID}}, acc::Set{GlobalID},
    visited::Set{GlobalID}, depth::Dict{GlobalID, Int}, 
    low::Dict{GlobalID, Int}, parent::Dict{GlobalID, GlobalID},
    vx::GlobalID, d::Int
)
    push!(visited, vx) 
    depth[vx] = d 
    low[vx] = d 
    num_children = 0 
    is_articulation = false

    for ni in edges[vx] 
        if ni in nodes 
            if !(ni in visited)
                parent[ni] = vx 
                articulation_points_ne!(nodes, edges, acc, visited, depth, low, parent, ni, d+1)
                num_children += 1
                if low[ni] >= depth[vx] 
                    is_articulation = true 
                end
                low[vx] = low[vx] < low[ni] ? low[vx] : low[ni]
            elseif !haskey(parent, vx) || ni != parent[vx] 
                low[vx] = low[vx] < depth[ni] ? low[vx] : depth[ni]
            end
        end
    end
    if (haskey(parent, vx) && is_articulation) || (!haskey(parent, vx) && num_children>1)
        push!(acc, vx) 
    end
end

function articulation_points(
    partition::MultiLevelPartition, 
    dists::Vector{Int}, 
    merged_nodes::Dict{GlobalID, Set{GlobalID}} = Dict{GlobalID, Set{GlobalID}}()
)::Set{GlobalID}
    acc = Set{GlobalID}()
    visited = Set{GlobalID}() 
    depth = Dict{GlobalID, Int}()
    low = Dict{GlobalID, Int}()
    parent = Dict{GlobalID, GlobalID}()
    vx_initial = any_pct_in(partition.subgraphs[dists[1]], partition.graph.root.global_ID)

    s = Stack{Tuple{GlobalID, Int}}()
    push!(s, (vx_initial, 0))
    #@show vx_initial 
    while length(s) > 0 
        vx, d = pop!(s) 

        if vx > 0 
            if vx in visited
                continue 
            end 
            push!(s, (-vx, -1))

            push!(visited, vx) 
            depth[vx] = d 
            low[vx] = d 
            num_children = 0 
          
            for (ni, MLN) in partition.graph.all_neighbors_map[vx] 
                if length(MLN.finer_neighbors) == 0 && look_up_district(partition, ni) in dists 
                    if vx in keys(merged_nodes) && ni in merged_nodes[vx] 
                        continue #skip the merged nodes adj to vx if vx is an art pt
                    end
                    if !(ni in visited)
                        parent[ni] = vx 
                        push!(s, (ni, d+1))
                        #num_children += 1 
                    #elseif vx == vx_initial || ni != parent[vx] 
                    #    low[vx] = low[vx] < depth[ni] ? low[vx] : depth[ni]  
                    end
                end
            end
            

        else 
            vx = -vx 
            # we've finished computing depth and low for all 
            # now we'll compute the lowpoint for vx and check whether vx is an articulation pt 

            # the low point of vx is the min of
            # - depth(vx) # this is the initial value of low[vx]
            # - depth(non-parent-neighbors(vx)) 
            # - low(children(vx))

            num_children = 0 
            for (ni, MLN) in partition.graph.all_neighbors_map[vx] 
                if length(MLN.finer_neighbors) == 0 && look_up_district(partition, ni) in dists 
                    if vx in keys(merged_nodes) && ni in merged_nodes[vx] 
                        continue #skip the merged nodes adj to vx if vx is an art pt
                    end

                    if ni !== vx_initial && parent[ni] == vx
                        num_children += 1
                        low[vx] = low[vx] < low[ni] ? low[vx] : low[ni]
                    end 
                    if vx !== vx_initial && parent[vx] !== ni 
                        low[vx] = low[vx] < depth[ni] ? low[vx] : depth[ni]  
                    end 

                    if ni !== vx_initial && parent[ni] == vx && # ni is a child of vx and
                        low[ni] >= depth[vx] && vx !== vx_initial # it satisfies the art-pt condition
                        push!(acc, vx)
                    end
                end
            end

            if vx == vx_initial && num_children > 1 
                push!(acc, vx) 
            end 
        end
    end 
    return acc 
end

function articulation_points_basic(
    partition::MultiLevelPartition, 
    dist::Int
)::Set{GlobalID}
    acc = Set{GlobalID}()
    visited = Set{GlobalID}() 
    depth = Dict{GlobalID, Int}()
    low = Dict{GlobalID, Int}()
    parent = Dict{GlobalID, GlobalID}()
    vx_initial = any_pct_in(partition.subgraphs[dist], partition.graph.root.global_ID)

    s = Stack{Tuple{GlobalID, Int}}()
    push!(s, (vx_initial, 0))
    while length(s) > 0 
        vx, d = pop!(s) 
        if vx > 0 
            if vx in visited continue end 
            push!(s, (-vx, -1))
            push!(visited, vx) 

            depth[vx] = d 
            low[vx] = d 
            num_children = 0 
          
            for ni in keys(partition.graph.all_neighbors_map[vx])
                if look_up_district(partition, ni) == dist
                    if !(ni in visited)
                        parent[ni] = vx 
                        push!(s, (ni, d+1))
                    end
                end
            end
        else 
            vx = -vx 
            num_children = 0 
            for ni in keys(partition.graph.all_neighbors_map[vx])
                if look_up_district(partition, ni) == dist
                    if ni !== vx_initial && parent[ni] == vx
                        num_children += 1
                        low[vx] = low[vx] < low[ni] ? low[vx] : low[ni]
                    end 
                    if vx !== vx_initial && parent[vx] !== ni 
                        low[vx] = low[vx] < depth[ni] ? low[vx] : depth[ni]  
                    end 

                    if ni !== vx_initial && parent[ni] == vx && # ni is a child of vx and
                        low[ni] >= depth[vx] && vx !== vx_initial # it satisfies the art-pt condition
                        push!(acc, vx)
                    end
                end
            end

            if vx == vx_initial && num_children > 1 
                push!(acc, vx) 
            end 
        end
    end 
    return acc 
end


# SNF will care about: d2n, n2d, subgraphs, 
#   dist_pops/areas/perims, adjacency, and adjacent_fine_neighbors
# SNF will clobber: shared_node, linking_edge_attrs 
function accept_SNF!(
    partition::MultiLevelPartition, 
    dist_to::Int, 
    pct::GlobalID; 
    merged_nodes = nothing, 
    known_artpts = Dict{Int, Set{GlobalID}}() # for perf and avoiding recomputation 
)
    dist_from = look_up_district(partition, pct)

    #d2n, n2d
    if merged_nodes !== nothing && pct in keys(merged_nodes) 
        for merged_node in merged_nodes[pct]
            mask_add_node!(partition.district_to_nodes[dist_to], partition.subgraphs[dist_to], merged_node) 
            mask_rem_node!(partition.district_to_nodes[dist_from], partition.subgraphs[dist_from], merged_node) 
            n2d_move_node!(partition, merged_node, dist_from, dist_to)    
        end
    else
        mask_add_node!(partition.district_to_nodes[dist_to], partition.subgraphs[dist_to], pct) 
        mask_rem_node!(partition.district_to_nodes[dist_from], partition.subgraphs[dist_from], pct) 
        n2d_move_node!(partition, pct, dist_from, dist_to)
    end

    #subgraphs 
    partition.subgraphs[dist_to] = take_subgraph(partition.graph, partition.district_to_nodes[dist_to])
    partition.subgraphs[dist_from] = take_subgraph(partition.graph, partition.district_to_nodes[dist_from])

    #pop, area 
    pct_pop = get_node_population(partition.graph, pct; merged_nodes=merged_nodes)
    pct_area = get_node_area(partition.graph, pct; merged_nodes=merged_nodes) 
    partition.dist_populations[dist_to] += pct_pop 
    partition.dist_populations[dist_from] -= pct_pop 
    partition.dist_areas[dist_to] += pct_area
    partition.dist_areas[dist_from] -= pct_area

    #perimeter 
    perim_with_dist_from = node_perim_with_district(partition, merged_nodes, pct, dist_from)
    perim_with_dist_to = node_perim_with_district(partition, merged_nodes, pct, dist_to)
    perim_total = node_perim(partition, pct; merged_nodes=merged_nodes)
    partition.dist_perimeters[dist_to] += (perim_total - 2*perim_with_dist_to)
    partition.dist_perimeters[dist_from] -= (perim_total - 2*perim_with_dist_from)

    # adjacency (and maybe shared nodes)
    #update_adjacency!(partition, dist_from, dist_to, pct)

    # adjacent_fine_neighbors
    if merged_nodes !== nothing && pct in keys(merged_nodes)
        for merged_node in merged_nodes[pct]
            # hopefully this works regardless of order. verify? 
            fix_adjacent_fine_neighbors!(partition, merged_node, dist_from, dist_to)
        end
    else
        fix_adjacent_fine_neighbors!(partition, pct, dist_from, dist_to)
    end

    # adjacency
    fix_adjacency!(partition)

    # articulation points; some perf changes to minimize recomputing art pts
    if dist_to in keys(known_artpts) && dist_from in keys(known_artpts)
        partition.articulation_points[dist_to] = known_artpts[dist_to]
        partition.articulation_points[dist_from] = known_artpts[dist_from]
    else
        partition.articulation_points[dist_to] = articulation_points_basic(partition, dist_to)
        partition.articulation_points[dist_from] = articulation_points_basic(partition, dist_from)
    end
end

function reject_SNF()
    return nothing
end


function dists_adj_to_node(partition, node_ID) 
    ds = Set()
    for (k,_) in partition.graph.all_neighbors_map[node_ID]
        push!(ds, look_up_district(partition, k))
    end
    delete!(ds, look_up_district(partition, node_ID))
    return ds 
end 

function adj_pcts(partition, node_ID) 
    fine_neighbors(partition.graph, node_ID)
    adj_nodes = [k for (k,_) in partition.graph.all_neighbors_map[node_ID]]
    return filter!(x -> length(get_children_IDs(partition.graph, x)) == 0, adj_nodes)
end

function fix_adjacent_fine_neighbors!(partition, node_ID, dist_from, dist_to) 
    adj_nodes = adj_pcts(partition, node_ID) 
    push!(adj_nodes, node_ID) # (me AND all my neighbors)
    adj_dists = Set([look_up_district(partition, x) for x in adj_nodes])
    
    # the only affected elements of adjacent_fine_neighbors are those where 
    # at least one of the indices is dist_to or dist_from. 
    # however i'll just be very coarse and remove all pairwise and add them all back 

    # remove all nodes from adj_fine_nbrs that live in the neighborhood of node_ID 
    # and whose connections involve dist_from or to
    for dist in adj_dists 
        if dist != -1 && dist !== dist_from 
            #@assert node_ID in partition.adjacent_fine_neighbors[dist, dist_from]
            delete!(partition.adjacent_fine_neighbors[dist, dist_from], node_ID)
        end
    end 
    for nbr in adj_nodes
        nbr_dist = look_up_district(partition, nbr)
        for dist in adj_dists
            if (dist == dist_from || nbr_dist == dist_from || 
                dist == dist_to || nbr_dist == dist_to) && (dist != -1 && nbr_dist != -1)
                delete!(partition.adjacent_fine_neighbors[dist, nbr_dist], nbr)
            end
        end
    end

    # now, for all the nodes <NBR> in the neighbordhood of node_ID, look at all the bordering 
    # districts of <NBR> and add <NBR> to those districts' adj_fine_nbrs 
    for nbr in adj_nodes
        nbr_dist = look_up_district(partition, nbr)
        dists_adj_to_nbr = dists_adj_to_node(partition, nbr) 
        for dist in dists_adj_to_nbr 
            if dist != -1 && nbr_dist != -1 && dist != nbr_dist
                push!(partition.adjacent_fine_neighbors[dist, nbr_dist], nbr)
            end 
        end
    end
end

# relies on adjacent_fine_neighbors being correct (i.e., updated)
function fix_adjacency!(partition) 
    for di=1:partition.num_dists 
        for dj = 1:partition.num_dists 
            partition.adjacency[di, dj] = length(partition.adjacent_fine_neighbors[di, dj]) > 0 
        end
    end
end

function node_perim_with_district(partition::MultiLevelPartition, 
    merged_nodes, node_ID::GlobalID, dist::Int; show=false)

    if merged_nodes !== nothing && node_ID in keys(merged_nodes) 
        perim = 0
        for merged_node in merged_nodes[node_ID]
            for (k,v) in partition.graph.all_neighbors_map[merged_node] 
                if !(k in merged_nodes[node_ID]) && 
                        length(v.finer_neighbors) == 0 && 
                        look_up_district(partition, k) == dist
                    perim += v.edge_attributes.length 
                end
            end
        end
        return perim 

    else
        perim = 0
        for (k,v) in partition.graph.all_neighbors_map[node_ID] 
            if length(v.finer_neighbors) == 0 && look_up_district(partition, k) == dist
                if show println("npwd", node_ID, ", ", dist, ", ", node_name(partition.graph, k)) end 
                perim += v.edge_attributes.length 
            end
        end
        return perim 
    end
end

function node_perim(partition::MultiLevelPartition, node_ID::GlobalID; 
    merged_nodes = nothing
) 
    perim = 0 
    for (k,v) in partition.graph.all_neighbors_map[node_ID] 
        if length(v.finer_neighbors) == 0 && 
                (merged_nodes === nothing || !(node_ID in keys(merged_nodes)) || !(k in merged_nodes[node_ID]))
            perim += v.edge_attributes.length 
        end
    end
    return perim 
end

# NB: this function does NOT work on multi level structures; only pct-county graphs. 
function n2d_move_node!(partition::MultiLevelPartition, node_ID::GlobalID, dist_from::Int, dist_to::Int) 
    n2d = partition.node_to_district

    @ensure !haskey(n2d, node_ID) || n2d[node_ID] == dist_from 
    node_parent = get_parent_ID(partition.graph, node_ID) 
    if !haskey(n2d, node_parent) || n2d[node_parent] == dist_from
        # node parent was intact. set parent county to -1  
        # and set all siblings to dist_from. 
        n2d[node_parent] = -1 
        foreach(x -> n2d[x] = dist_from, get_children_IDs(partition.graph, node_parent))
        n2d[node_ID] = dist_to 
    elseif n2d[node_parent] == -1 
        # it was split before. need to check whether it's still split. 
        n2d[node_ID] = dist_to 
        sibling_dists = Set([n2d[x] for x in get_children_IDs(partition.graph, node_parent)])
        if length(sibling_dists) > 1 || -1 in sibling_dists 
            # it is still split (and was split before). do nothing 
        else 
            # it is now combined. all children of node_parent belong to dist_to now. 
            # n2d[node_parent] = dist_to; remove all siblings of node_ID from n2d 
            @assert dist_to != -1
            n2d[node_parent] = dist_to 
            foreach(x -> delete!(n2d, x), get_children_IDs(partition.graph, node_parent))
        end
    else 
        @assert false "n2d not consistent with SNF pct"
    end
end

# adds the given node to the mask
# no guarantees if you try to add a node already in the mask  
function mask_add_node!(mask::Mask, graph::MultiLevelGraph, node_ID::GlobalID; with=nothing) 
    #println("mask add: ", node_ID) 
    @assert mask == graph.mask 
    if is_in_graph(graph, node_ID)
        @show node_ID 
    end 
    @assert !(is_in_graph(graph, node_ID))
    node_parent = get_parent_ID(graph, node_ID) 
    if haskey(mask, node_parent) 
        @assert mask[node_parent] !== nothing 
    else 
        mask_add_node!(mask, graph, node_parent; with=Set()) 
    end
    push!(mask[node_parent], node_ID)
    mask[node_ID] = with
  
    # all of node_parent's children are in the mask...
    my_siblings = get_children_IDs(graph.parent, node_parent)
    if length(mask[node_parent]) == length(my_siblings)
        if all([mask[x] === nothing for x in my_siblings])
        #.. and all of them are fully intact. sparse rep means := `nothing` 
            #println("county fully un-split")
            mask[node_parent] = nothing 
            for x in my_siblings
                delete!(mask, x)
            end
        end 
    end
end

function mask_expand_ancestors!(mask::Mask, graph::MultiLevelGraph, node::GlobalID) 
    # nb: ancestors includes self 
    if node !== graph.root.global_ID && (!haskey(mask, node) || mask[node] === nothing)
        mask[node] = get_children_IDs(graph, node) 
        mask_expand_ancestors!(mask, graph, get_parent_ID(graph, node))
    end 
end

# removes the given node from the mask. if it does not exist, errors. 
# no guarantees if you try to remove the root or the last node in the mask 
# mask has to correspond to the given graph (obviously, hopefully) 
# this is incorrect; graph should be the base graph that the mask subgraphs from 
function mask_rem_node!(mask::Mask, graph::MultiLevelGraph, node_ID::GlobalID)
    @assert is_in_graph(graph, node_ID)
    #println("mask remove: ", node_ID) 
    node_parent = get_parent_ID(graph, node_ID) 
    if !haskey(mask, node_parent)

        mask_expand_ancestors!(mask, graph.parent, node_parent)
    elseif mask[node_parent] === nothing 
        #println("splitting new county")
        mask[node_parent] = Set(get_children_IDs(graph.parent, node_parent))
        for child in mask[node_parent] 
            mask[child] = nothing
        end
    end
    pop!(mask[node_parent], node_ID)
    delete!(mask, node_ID) 

    if length(mask[node_parent]) == 0 
        if node_parent == graph.root.global_ID
            @assert false "you seem to be deleting the last node in the district. this might be because your pop deviations allow 0-population districts."
        end
        mask_rem_node!(mask, graph, node_parent)
    end
end

function build_single_node_flip(constraints)
    function f(partition, measure, rng, i) 
        
        p, reject, accept, pct, dist_to = single_node_flip_basic(partition, measure, constraints, rng)

        p_back = single_node_flip_basic(partition, measure, constraints, rng, 
            just_prob_of_flip = (dist_to, pct)
        )
    end
end

function single_node_flip_basic(
    partition::MultiLevelPartition,
    measure::Measure, # this is the measure we'll TEMPER on. 
    constraints::Dict, #todo: deal with constraints? 
    rng::AbstractRNG;
    just_prob_of_flip::Union{Nothing, Tuple{Int, GlobalID}} = nothing, 
    extra::Any = nothing 
)
    alpha = 1.0 # tempering parameter 

    merged_nodes = Dict{GlobalID, Set{GlobalID}}()
    
    #compute the articulation points (which cannot be flipped lest the district they're in become disconnected)
    art_pts = Dict{Int, Set{GlobalID}}()
    for d = 1:partition.num_dists  
        art_pts[d] = partition.articulation_points[d] 
        # what???
        #articulation_points(partition, [d], merged_nodes)
    end
    

    # compute all potential flips 
    potential_flips = Dict{Int, Set{GlobalID}}() 
    for district_from = 1:partition.num_dists 
        for district_to = 1:partition.num_dists     
            adjacent_pcts = partition.adjacent_fine_neighbors[district_to, district_from]
            

            function is_valid(x) 
                if x in art_pts[district_from] 
                    return false 
                end 
                xpop = get_node_population(partition.graph, x) 
                return partition.dist_populations[district_to] + xpop <= constraints[PopulationConstraint].max_pop && 
                    partition.dist_populations[district_from] - xpop >= constraints[PopulationConstraint].min_pop
            end 

            this_valid_flips = Set(filter(is_valid, collect(adjacent_pcts)))
            if haskey(potential_flips, district_to) 
                union!(potential_flips[district_to], this_valid_flips)
            else 
                potential_flips[district_to] = this_valid_flips
            end
        end
    end

    #@show sum([length(x) for (_,x) in potential_flips]) 

    initial_energies =partially_compute_energy(partition, measure)
    og_energy = fully_compute_energy(partition, measure)

    # compute the energy for each potential node flip 
    potential_energies = Dict{Int, Dict{GlobalID, Real}}([k => Dict{GlobalID, Real}() for k in keys(potential_flips)])
    total_energies = 0
    npcs = 0

    #og_iso = sum((partition.dist_perimeters .^2) ./ partition.dist_areas)
    for (dist_to, pcts) in potential_flips 
        for pct in pcts 
            npcs += 1 
            this_energy = finish_energy_with_SNF(initial_energies, partition, merged_nodes, measure, dist_to, pct)
            #this_energy = finish_energy_with_SNF(partition, measure.weights[1], dist_to, pct, og_iso)
            potential_energies[dist_to][pct] = exp(-alpha*(this_energy - og_energy))
			if isnan(potential_energies[dist_to][pct]) 
				potential_energies[dist_to][pct] = 0 
			end
            total_energies += potential_energies[dist_to][pct]
        end
    end
    
    @assert npcs > 0 "there are no valid single node flips?"

    if just_prob_of_flip !== nothing
        dist_to, pct_to_flip = just_prob_of_flip 
        @assert haskey(potential_energies, dist_to)
        if !haskey(potential_energies[dist_to], pct_to_flip)
            println("--- would be asserting false here... ---")
            println(dist_to, " ", pct_to_flip, " ", length(partition.graph.ID_to_node))
            #return 0
            @assert false 
        end
        res = potential_energies[dist_to][pct_to_flip] / total_energies 
        return res

    else 
        #write_energies_to_file(partition, potential_energies, "viz/energy.csv") 

        # choose proportional to the energy 
        choice = rand(rng) * total_energies 
        running_total = 0 
        for (dist_to, pcts) in potential_flips
            for pct in pcts 
                running_total += potential_energies[dist_to][pct]
                if running_total > choice 
                    #@show pct, dist_to, dist_from
                    reject_update = (_ -> nothing)
                    accept_update = (new_aps -> accept_SNF!(partition, dist_to, pct; merged_nodes=merged_nodes, known_artpts = new_aps))

                    p_fwd = potential_energies[dist_to][pct] / total_energies

                    return p_fwd, reject_update, accept_update, -1/alpha * log(potential_energies[dist_to][pct]), pct, dist_to
                end
            end
        end
        @show total_energies 
        @show potential_energies 
        @show fully_compute_energy(partition, measure)
        @show extra 
        @assert false 
    end
end
