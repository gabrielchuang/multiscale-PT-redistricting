# HierarchicalPartition: Partition structure for perfect-matching-based districtings. 

struct HierarchicalPartition 
    graph::MultiLevelGraph # just the precinct graph. 
        #Question: how to think about county splitting? Can incorporate into the energy? 
   
    oriented_neighbors_by_level::Vector{Dict{GlobalID, Vector{GlobalID}}}
        # this doesn't really have to be *by level*; it could just be Dict{GlobalID, Vector{GlobalID}} too

    node_parents::Dict{GlobalID, GlobalID}
    node_children::Dict{GlobalID, Set{GlobalID}} # or maybe Union{GlobalID, Tuple{GlobalID}}?  

    nodes_by_level::Vector{Set{GlobalID}}
    edges_by_level::Vector{Dict{GlobalID, Set{GlobalID}}}

    node_attributes::Dict{GlobalID, Dict{String, Any}}
    edge_attributes::Dict{UndirectedEdge, EdgeAttributes}
end

#= 
Things we need to be able to do: 
- adjacency matrix + oriented neighbors at a given level (e.g., to construct H4 from H3)
- adjacency matrix + oriented neighbors for some subset, at a given level 
    (e.g., to reconstruct H3 given H4 and H2 means getting the size-4 subgraphs of each H4)
- compute area/perim/pop/etc of coarser nodes given the finer nodes they're made of  

oriented_nbrs = JSON.parsefile("test_graphs/pct21_oriented_nbrs.json")
oriented_nbrs = Dict{GlobalID, Vector{GlobalID}}([parse(Int, k) => v for (k,v) in oriented_nbrs])
=# 
 
function HierarchicalPartition(G::MultiLevelGraph, onbrs::Dict{GlobalID, Vector{GlobalID}})
    nodes = finest_level_nodes(G)
    edges = finest_level_edges_adjlist(G)
    edge_attrs = Dict(
        [UndirectedEdge(x,y) => G.all_neighbors_map[x][y].edge_attributes 
            for x in nodes for y in edges[x]])

    for node in nodes
        node_attributes(G, node)
    end

    return HierarchicalPartition(
        G, 
        [onbrs], 
        Dict(), # no node parents
        Dict(), # no node children 
        [nodes], 
        [edges],
        G.node_attributes, 
        edge_attrs
    )
end

function topmost_ancestor(hp::HierarchicalPartition, node_ID::GlobalID) 
    if haskey(hp.node_parents, node_ID) 
        return topmost_ancestor(hp, hp.node_parents[node_ID])
    else
        return node_ID 
    end
end

function add_level_at_top!(
    hp::HierarchicalPartition, 
    rng, 
    matching;
    unmatched_nodes = [],
    exclude_nodes = nothing
)
    #nodes = hp.nodes_by_level[end]
    
    #if exclude_nodes !== nothing
    #    nodes = filter(x -> !exclude_nodes(x), nodes)
    #end
    #@assert mod(length(nodes), 2) == 0
    #QEM = QuadEdgeMesh(nodes, hp.oriented_neighbors_by_level[end])
    #matching = fast_sample_pm(nodes, hp.edges_by_level[end], QEM, rng, (_,_)->1)

    max_node_ID = maximum(maximum.(hp.nodes_by_level))
    new_node_ID = max_node_ID + 1

    next_level_onbrs = Dict{GlobalID, Vector{Vector{GlobalID}}}()
    # maps a node to the set of runs of onbrs. 

    new_level_num = length(hp.nodes_by_level) + 1
    push!(hp.nodes_by_level, Set{GlobalID}())

    for matched_pair in matching 
        # set: 
        # node parents/children
        # oriented nbrs 
        # node and edge attrs 
        c1, c2 = matched_pair.id1, matched_pair.id2
        hp.node_parents[c1] = new_node_ID
        hp.node_parents[c2] = new_node_ID
        hp.node_children[new_node_ID] = Set{GlobalID}([c1,c2])

        # combine the two nodes
        push!(hp.nodes_by_level[new_level_num], new_node_ID)
        hp.node_attributes[new_node_ID] =  merge_node_attrs([hp.node_attributes[c1], hp.node_attributes[c2]])

        next_level_onbrs[new_node_ID] = combine_onbrs(c1, c2, hp.oriented_neighbors_by_level[end])
        new_node_ID += 1 
    end

    for unmatched_node in unmatched_nodes 
        hp.node_parents[unmatched_node] = new_node_ID 
        hp.node_children[new_node_ID] = Set{GlobalID}(unmatched_node)
        push!(hp.nodes_by_level[new_level_num], new_node_ID)
        hp.node_attributes[new_node_ID] =  hp.node_attributes[unmatched_node]
        next_level_onbrs[new_node_ID] = [hp.oriented_neighbors_by_level[end][unmatched_node]]
        new_node_ID += 1
    end

    coarsened_onbrs = Dict{GlobalID, Vector{GlobalID}}()
    for (k,v) in next_level_onbrs
        coarse_onbrs = []
        for section in v
            this_section_parents = [x == -1 ? -1 : hp.node_parents[x] for x in section]
            remove_adj_duplicates_cyclic!(this_section_parents)
            push!(coarse_onbrs, this_section_parents)
        end
        coarsened_onbrs[k] = collect(flatten(coarse_onbrs))
    end
    push!(hp.oriented_neighbors_by_level, coarsened_onbrs)

    # need to fill: edges by level and edge attributes 
    edges_at_new_level = Dict{GlobalID, Set{GlobalID}}()
    for (x1, nbrs) in hp.edges_by_level[end]
        p1 = hp.node_parents[x1] 
        edges_at_new_level[p1] = Set()
    end
    for (x1, nbrs) in hp.edges_by_level[end]
        p1 = hp.node_parents[x1] 
        for x2 in nbrs 
            if x1 < x2 continue end #this prevents us from adding every edge twice 
            p2 = hp.node_parents[x2]
            if p1 == p2 continue end 
            push!(edges_at_new_level[p1], p2)
            push!(edges_at_new_level[p2], p1)
            p_edge = UndirectedEdge(p1, p2)
            if haskey(hp.edge_attributes, p_edge)
                hp.edge_attributes[p_edge] = combine_edge_attrs(
                    hp.edge_attributes[p_edge], 
                    hp.edge_attributes[UndirectedEdge(x1, x2)])
                    # note: this adds every edge twice!! 
            else 
                hp.edge_attributes[p_edge] = hp.edge_attributes[UndirectedEdge(x1, x2)]
            end
        end
    end

    push!(hp.edges_by_level, edges_at_new_level)
end

function remove_level_at_top!(hp::HierarchicalPartition)
    for node in hp.nodes_by_level[end] 
        for child in hp.node_children[node]
            delete!(hp.node_parents, child)
        end
        if haskey(hp.edges_by_level[end], node)
            for n2 in hp.edges_by_level[end][node]
                delete!(hp.edge_attributes, UndirectedEdge(node, n2))
            end
        end
        delete!(hp.node_children, node)
        delete!(hp.node_attributes, node)
    end
    deleteat!(hp.edges_by_level, length(hp.edges_by_level))
    deleteat!(hp.nodes_by_level, length(hp.nodes_by_level))
    deleteat!(hp.oriented_neighbors_by_level, length(hp.oriented_neighbors_by_level))
end

function splice_at(onbrs, at)
    i = findfirst(x -> x == at, onbrs) 
    return vcat(onbrs[i+1:end], onbrs[1:i-1])
end

function multisplice_at(onbrs, c2)
    sections = [[]] 
    for nbr in onbrs 
        if nbr == c2 
            push!(sections, []) 
        else 
            push!(sections[end], nbr)
        end
    end
    sections[1] = vcat(sections[end], sections[1])
    deleteat!(sections, length(sections))
    return sections
end

function remove_adj_duplicates_cyclic!(v)
    if length(v) <= 1 return end 
    x = v[1] 
    i = 2 
    while i <= length(v) 
        if v[i] == x 
            deleteat!(v, i) 
        else 
            x = v[i]
            i += 1 
        end
    end
    if v[1] == v[end] && length(v) > 1
        deleteat!(v, length(v))
    end
end


# kind of a gross function with lots of cases. 
# Given the oriented neighbors, compute the oriented neighbors of the node formed if 
# you merged them. 
function combine_onbrs(
    c1::GlobalID,
    c2::GlobalID,
    all_onbrs::Dict{GlobalID, Vector{GlobalID}}
)
    onbrs1 = all_onbrs[c1]
    onbrs2 = all_onbrs[c2]
    if count(x-> x == c2, onbrs1) == 1 && count(x-> x == c1, onbrs2) == 1
        return [vcat(splice_at(onbrs1, c2), splice_at(onbrs2, c1))]
    else 
        if count(x-> x == c2, onbrs1) != count(x-> x == c1, onbrs2)
            @show c1, c2 
            @assert false 
        end 
        # get the sections of nbrs of c1 split by c2 
        sections1 = multisplice_at(onbrs1, c2) 
        sections2 = multisplice_at(onbrs2, c1) 
        @assert length(sections1) == length(sections2)

        onbrs = []
        # pair up the sections 
        for (i1, s1) in enumerate(sections1)
            succ = false 
            for (i2, s2) in enumerate(sections2)
                if length(intersect(Set(s1), Set(s2))) > 0 
                    # ok, nice, these two runs of onbrs intersect. They must 
                    # merge into one run, then. 
                    push!(onbrs, vcat(s1, s2))
                    deleteat!(sections2, i2)
                    succ = true 
                end
            end
            if !succ 
                # uh oh 
                println("a slightly concerning case... (onbr merging)")
                s1_end_nbrs = union(Set(all_onbrs[s1[1]]), Set(all_onbrs[s1[end]]))
                for (i2, s2) in enumerate(sections2)
                    if length(intersect(Set(s1_end_nbrs), Set(s2))) > 0 
                        # ok, one of the neighbors of the ends of s1 is an onbr of s2. 
                        push!(onbrs, vcat(s1, s2))
                        deleteat!(sections2, i2)
                        succ = true 
                    end
                end
            end
            if !succ 
                # oh well... 
                # somehow ths run of onbrs doesn't touch anything in s2. I don't really 
                # know how this could happen. 
                println("isolated oriented neighbor that somehow borders both sides??")
                push!(onbrs, s1) 
            end 
        end
        for remaining_s2 in sections2 
            # same as the case above but for s2 
            println("isolated oriented neighbor that somehow borders both sides??")     
            push!(onbrs, remaining_s2)
        end
        for section in onbrs 
            remove_adj_duplicates_cyclic!(section)
        end
        #@show sections1, sections2 
        #@show c1, c2 
        #@show onbrs 
        return onbrs 
    end
end

function check_onbr_numbers_match(onbrs) 
    for (k,v) in onbrs
        for y in v
            if y == -1 continue end 
            if count(x-> x == y, v) != count(x-> x == k, onbrs[y])
                onbrsy = onbrs[y]
                println("WARNING: $k and $y have mismatched onbrs: $v, $onbrsy")
                if count(x-> x == y, v) > count(x-> x == k, onbrs[y])
                    deleteat!(v, findfirst(x->x==y, v))
                else
                    @assert count(x-> x == y, v) < count(x-> x == k, onbrs[y])
                    deleteat!(onbrs[y], findfirst(x->x==k, onbrs[y]))
                end 
                println("NEW onbrs (hacky; may be wrong): $v, $onbrsy")
            end
        end
    end
end

function node_pop(hp::HierarchicalPartition, node_ID) 
    return hp.node_attributes[node_ID][hp.graph.base_graph.pop_col]
end

using RandomNumbers


function add_level_with_pop_threshold!(
    hp::HierarchicalPartition, 
    pop_threshold, 
    rng
)
    nodes = hp.nodes_by_level[end]
    edges = deepcopy(hp.edges_by_level[end])

    nodes = filter(x -> node_pop(hp, x) <= pop_threshold, nodes)
    println("# nodes below pop threshold: ", length(nodes))
    filter_edges!(edges, nodes)
    ccs = even_ccs(nodes, edges) 
    matching = []
    unmatched_nodes = deepcopy(hp.nodes_by_level[end])
    log_p_total = 0

    for cc_nodes in ccs
        if length(cc_nodes) == 0 continue end
        cc_edges = deepcopy(edges)
        filter_edges!(cc_edges, cc_nodes)

        print(length(cc_nodes), " ")
        
        if length(cc_nodes) == 910 
            @show cc_nodes 
        end
        
        QEM0 = QuadEdgeMesh(cc_nodes, hp.oriented_neighbors_by_level[end])
        remove_multiedges!(QEM0)
        #@show QEM0.edge_lookup 

        rng2 = PCG.PCGStateOneseq(UInt64, 15152)
        matching0, log_p = fast_sample_pm(cc_nodes, cc_edges, QEM0, rng, 
            (src, dst) -> #hp.edge_attributes[UndirectedEdge(src, dst)].length
               maximum([
                10*(1 - abs(node_pop(hp, src) + node_pop(hp, dst) - pop_threshold) / pop_threshold)^2, 
                0.01])
            )
        log_p_total += log_p 
        append!(matching, matching0)
        if length(matching0) !== 0
            setdiff!(unmatched_nodes, cc_nodes)
        end
    end

    println("# nodes that got paired: ", length(hp.nodes_by_level[end]) - length(unmatched_nodes))

    add_level_at_top!(hp, rng, matching; unmatched_nodes = unmatched_nodes)
    return log_p_total 
end

function add_level_with_pop_target!(
    hp::HierarchicalPartition, 
    pop_target, 
    rng
)
    nodes = hp.nodes_by_level[end]
    @assert mod(length(nodes), 2) == 0 
    edges = deepcopy(hp.edges_by_level[end])

    println("# nodes to match: ", length(nodes))

    QEM0 = QuadEdgeMesh(nodes, hp.oriented_neighbors_by_level[end])
    remove_multiedges!(QEM0)

    matching, log_p = fast_sample_pm(nodes, edges, QEM0, rng, 
        (src, dst) -> weight(node_pop(hp, src) + node_pop(hp, dst), pop_target, 1))

    add_level_at_top!(hp, rng, matching; unmatched_nodes = [])
    return log_p
end

function weight(p, p_target, beta)
    if p > p_target 
        return exp(-beta*(p/p_target)^2)
    else 
        return exp(-beta*(p_target/p)^2)
    end
end

function pm_count(hp::HierarchicalPartition, rng) 
    log_pmc = 0
    for di in hp.nodes_by_level[end]
        nodes = Set{GlobalID}()
        add_descendants_at_level!(nodes, hp, di, length(hp.nodes_by_level), 1)

        vmap = reverse_vmap(collect(nodes))

        edges = deepcopy(hp.edges_by_level[1])
        filter_edges!(edges, nodes)

        QEM = QuadEdgeMesh(nodes, hp.oriented_neighbors_by_level[1])
        remove_multiedges!(QEM)
    
        edge_set = Set{UndirectedEdge}([UndirectedEdge(x,y) for (x, v) in edges for y in v])
        vmap, tutte::Matrix{Float64} = construct_pfaffian_orientation(edge_set, nodes, QEM, rng, (_, _) -> 1)

        log_pmc += 1/2 * logdet(tutte)
    end
    return log_pmc
end

#this one is approximate bc a non hierarhcical MLP can have odd # nodes or no PMS etc 
function pm_count(partition::MultiLevelPartition, onbrs, rng) 
    log_pmc = 0
    nskips = 0
    for di = 1:partition.num_dists
        try 
            nodes = deepcopy(partition.district_to_nodes[di][partition.graph.root.global_ID])
            if mod(length(nodes), 2) == 1
                delete!(nodes, first(nodes))
            end
            vmap = reverse_vmap(collect(nodes))

            g = partition.graph 
            edges = Dict([k => collect(keys(v)) for (k,v) in g.all_neighbors_map])
            filter_edges!(edges, nodes)

            QEM = QuadEdgeMesh(nodes, onbrs)
            remove_multiedges!(QEM)
        
            edge_set = Set{UndirectedEdge}([UndirectedEdge(x,y) for (x, v) in edges for y in v])
            vmap, tutte::Matrix{Float64} = construct_pfaffian_orientation(edge_set, nodes, QEM, rng, (_, _) -> 1)

            factor = 1/2 * logdet(tutte)
            if factor == -Inf 
                nskips += 1 
                continue 
            end
            log_pmc += factor
        catch e 
            @show e 
        end
    end
    log_pmc += nskips * (log_pmc / (partition.num_dists - nskips))
    return log_pmc
end

function spanning_forest_count(hp::HierarchicalPartition) 
    log_sfc = 0
    for di in hp.nodes_by_level[end]
        nodes = Set{GlobalID}()
        add_descendants_at_level!(nodes, hp, di, length(hp.nodes_by_level), 1)
        vmap = reverse_vmap(collect(nodes))

        edges = deepcopy(hp.edges_by_level[1])
        filter_edges!(edges, nodes)

        srcs = Vector{Int64}() 
        dsts = Vector{Int64}()
        for (x1, nbrs) in edges 
            for x2 in nbrs 
                push!(srcs, vmap[x1]) 
                push!(dsts, vmap[x2]) 
            end
        end

        g = SimpleWeightedGraph(srcs, dsts, ones(length(srcs)))
        log_sfc += log_nspanning(g)
    end
    return log_sfc 
end

# abs(node_pop(hp, src) + node_pop(hp, dst) - pop_target)

function isoperimetric_score(hp, n1, n2, level) 
    l1 = sum([hp.edge_attributes[UndirectedEdge(n1, x)].length for x in hp.edges_by_level[level][n1]])
    l2 = sum([hp.edge_attributes[UndirectedEdge(n2, x)].length for x in hp.edges_by_level[level][n2]])
    l12 = hp.edge_attributes[UndirectedEdge(n1,n2)].length
    a1 = hp.node_attributes[n1]["area"]
    a2 = hp.node_attributes[n2]["area"]

    return (l1 + l2 - l12 - l12)^2 / (a1 + a2)
end

function add_descendants_at_level!(descs, hp, x, ancestor_level, target_level)
    @assert x in hp.nodes_by_level[ancestor_level]
    if ancestor_level == target_level 
        push!(descs, x)
    else 
        for y in hp.node_children[x]
            add_descendants_at_level!(descs, hp, y, ancestor_level-1, target_level)
        end
    end 
end

function g_ancestor(hp, node_ID, level=1)
    if level == length(hp.nodes_by_level)
        @assert node_ID in hp.nodes_by_level[end]
        return node_ID
    else 
        return g_ancestor(hp, hp.node_parents[node_ID], level+1)
    end
end


function add_level_with_edge_reweight!(
    hp::HierarchicalPartition, 
    tuttes_by_level, 
    rng, 
    beta=1, # a factor to change how much we modify the weights to account for Z_PM 
    iso_weight=0.1 # an iso weight 
)
    @assert length(tuttes_by_level) == length(hp.nodes_by_level) - 1
    curr_level = length(hp.nodes_by_level)

    function this_level_edge_weights(n1, n2) 
        weight = exp(-iso_weight * isoperimetric_score(hp, n1, n2, curr_level))
        if beta != 0
            for l = 1:curr_level-1
                these_nodes = Set{GlobalID}() 
                add_descendants_at_level!(these_nodes, hp, n1, curr_level, l)
                add_descendants_at_level!(these_nodes, hp, n2, curr_level, l)
                vmap, tutte = tuttes_by_level[l]
                vmapped_nodes = [vmap[x] for x in these_nodes]
                PM_count = sqrt(det(tutte[vmapped_nodes, vmapped_nodes]))
                weight /= exp(beta * log(PM_count)) 
            end
        end
        return weight
    end

    # w_i = e^-J / prod l=1 to i-1 of ZPM(...)

    nodes = deepcopy(hp.nodes_by_level[end])
    QEM = QuadEdgeMesh(nodes, hp.oriented_neighbors_by_level[end])
    remove_multiedges!(QEM)

    edges = deepcopy(hp.edges_by_level[end])
    @assert all([length(v) !== 0 for (k,v) in edges])

    edge_set = Set{UndirectedEdge}([UndirectedEdge(x,y) for (x, v) in edges for y in v])
    vmap, tutte::Matrix{Float64} = construct_pfaffian_orientation(edge_set, nodes, QEM, rng, this_level_edge_weights)
    push!(tuttes_by_level, (vmap, deepcopy(tutte)))
    @assert all([length(v) !== 0 for (k,v) in edges])

    matching, log_p = fast_sample_pm(nodes, edges, vmap, tutte, rng)
    @assert length(matching) == length(hp.nodes_by_level[end]) / 2
    add_level_at_top!(hp, rng, matching; unmatched_nodes = [])
    return log_p 
end

# just makes a matching with iso weight for the given graph (and its oriented nbrs...)
function iso_matching(Gnc, oriented_nbrs, target_pop, pop_weight, rng) 
    hp = HierarchicalPartition(Gnc, oriented_nbrs) 
    nodes = deepcopy(hp.nodes_by_level[end])
    QEM = QuadEdgeMesh(nodes, hp.oriented_neighbors_by_level[end])

    remove_multiedges!(QEM)
    edges = deepcopy(hp.edges_by_level[end])
    edge_set = Set{UndirectedEdge}([UndirectedEdge(x,y) for (x, v) in edges for y in v])

    function edge_weight_function(n1, n2) 
        iso_weight = 0.02
        weight = exp(-iso_weight * isoperimetric_score(hp, n1, n2, 1))
        
        npop = get_node_population(Gnc, n1) + get_node_population(Gnc, n2)
        dist_from_target = npop < target_pop ? -5 : 3*(npop/target_pop)
        weight *= exp(-pop_weight * dist_from_target)
        
        return weight
    end

    vmap, tutte::Matrix{Float64} = construct_pfaffian_orientation(edge_set, nodes, QEM, rng, edge_weight_function)
    matching, _ = fast_sample_pm(nodes, edges, vmap, tutte, rng)
    return matching 
end

function iso_matching3(Gnc, oriented_nbrs, target_pop, pop_weight, rng) 
    hp = HierarchicalPartition(Gnc, oriented_nbrs) 
    nodes = deepcopy(hp.nodes_by_level[end])
    edges = deepcopy(hp.edges_by_level[end])
    edge_set = Set{UndirectedEdge}([UndirectedEdge(x,y) for (x, v) in edges for y in v])
    
    if mod(length(nodes), 2) == 1
        delnode = collect(nodes)[5]
        @show node_name(Gnc, delnode), delnode
        delete!(nodes, delnode)
        for x in edge_set 
            if x.id1 == delnode || x.id2 == delnode 
                delete!(edge_set, x)
                delete!(edges[x.id2], x.id1)
                delete!(edges[x.id1], x.id2)
            end
        end
        delete!(edges, delnode)
    end

    QEM = QuadEdgeMesh(nodes, hp.oriented_neighbors_by_level[end])
    remove_multiedges!(QEM)



    function edge_weight_function(n1, n2) 
        iso_weight = 0.2
        weight = exp(-iso_weight * isoperimetric_score(hp, n1, n2, 1))
        
        npop = get_node_population(Gnc, n1) + get_node_population(Gnc, n2)
        dist_from_target = npop < target_pop ? -5 : 3*(npop/target_pop)
        weight *= exp(-pop_weight * dist_from_target)

        if ((0.5 * get_node_perimeter(Gnc, n1)^2 / get_node_area(Gnc, n1) + 
            0.5 * get_node_perimeter(Gnc, n2)^2 / get_node_area(Gnc, n2) - 
            isoperimetric_score(hp, n1, n2, 1)) < -15)
            weight = exp(-20)
        end

        
        return weight
    end

    vmap, tutte::Matrix{Float64} = construct_pfaffian_orientation(edge_set, nodes, QEM, rng, edge_weight_function)

    matching, _ = fast_sample_pm(nodes, edges, vmap, tutte, rng)
    return matching, edge_weight_function 
end

function iso_matching2(Gnc, oriented_nbrs, target_pop, pop_weight, rng) 
    hp = HierarchicalPartition(Gnc, oriented_nbrs) 
    nodes = deepcopy(hp.nodes_by_level[end])
    QEM = QuadEdgeMesh(nodes, hp.oriented_neighbors_by_level[end])
    remove_multiedges!(QEM)
    edges = deepcopy(hp.edges_by_level[end])
    edge_set = Set{UndirectedEdge}([UndirectedEdge(x,y) for (x, v) in edges for y in v])

    function edge_weight_function(n1, n2) 
        iso_weight = 0.1
        weight = exp(-iso_weight * isoperimetric_score(hp, n1, n2, 1))
        
        npop = get_node_population(Gnc, n1) + get_node_population(Gnc, n2)
        dist_from_target = (abs(npop - target_pop) / target_pop) 
            #npop < target_pop ? (target_pop/npop) : (npop/target_pop)
        weight *= exp(-pop_weight * dist_from_target)

        if count(x -> x == -1, combine_onbrs(n1, n2, oriented_nbrs)) > 1 
            weight = exp(-200)
        end
        
        return weight
    end

    vmap, tutte::Matrix{Float64} = construct_pfaffian_orientation(edge_set, nodes, QEM, rng, edge_weight_function)
    matching, _ = fast_sample_pm(nodes, edges, vmap, tutte, rng)
    return matching 
end

function matching_maximize_below_thres(Gnc, oriented_nbrs, thres, rng) 
    hp = HierarchicalPartition(Gnc, oriented_nbrs) 
    nodes = deepcopy(hp.nodes_by_level[end])
    QEM = QuadEdgeMesh(nodes, hp.oriented_neighbors_by_level[end])
    remove_multiedges!(QEM)
    edges = deepcopy(hp.edges_by_level[end])
    edge_set = Set{UndirectedEdge}([UndirectedEdge(x,y) for (x, v) in edges for y in v])

    function edge_weight_function(n1, n2) 
        npop = get_node_population(Gnc, n1) + get_node_population(Gnc, n2)
        weight = npop < thres ? 4 : 1
        
        return weight
    end

    vmap, tutte::Matrix{Float64} = construct_pfaffian_orientation(edge_set, nodes, QEM, rng, edge_weight_function)
    matching, _ = fast_sample_pm(nodes, edges, vmap, tutte, rng)
    return matching 
end