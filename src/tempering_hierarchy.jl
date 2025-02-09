# a hierarchy for holding information for swapping up and down 
# and some functions for constructing and accessing such hierarchies 

# invariant: 
# all parents with only one child have the same node_ID as their child. 
struct TemperingHierarchy 
    graphs_by_level::Vector{MultiLevelGraph}
    #graph_files_by_level::Vector{String}
    #matching::Vector{UndirectedEdge} # i guess... idk. but just for now 
        # we'll just have this hierarchy only deal with the levels from 
        # base to pairs, or similar. 

    # |children_by_level| = |graphs_by_level| 
    # |parents_by_level| = |graphs_by_level| - 1
    children_by_level::Vector{Dict{GlobalID, Vector{GlobalID}}}
    parents_by_level::Vector{Dict{GlobalID, GlobalID}}
    onbrs_by_level::Vector{Dict{GlobalID, Vector{GlobalID}}}
    filename_prefix::String
end

function TemperingHierarchy(G, onbrs, prefix::String) 
    return TemperingHierarchy([G], [Dict{GlobalID, Vector{GlobalID}}()], [], [onbrs], prefix) 
end

function add_level!(TH, matching)
    new_level_num = length(TH.graphs_by_level)+1

    # do the pairing up 
    new_filename = TH.filename_prefix 
    pair_parent, pair_children, new_onbrs = 
        pair_up_nodes(TH.graphs_by_level[end], matching, new_filename, new_level_num,
        TH.onbrs_by_level[end])
    push!(TH.parents_by_level, pair_parent)
    push!(TH.children_by_level, pair_children) 
    push!(TH.onbrs_by_level, new_onbrs)

    # make the new graph 
    #node_data = Set(["block_name", "pop2020cen", "area", "border_length", "id"])
    #new_base_graph = BaseGraph(new_filename*string(new_level_num)*".json", "pop2020cen", inc_node_data=node_data, area_col="area",node_border_col="border_length", edge_perimeter_col="length")
    node_data = Set(["block_name", TH.graphs_by_level[1].base_graph.pop_col, "area", "border_length", "id"])
    new_base_graph = BaseGraph(new_filename*string(new_level_num)*".json", TH.graphs_by_level[1].base_graph.pop_col, inc_node_data=node_data, area_col="area",node_border_col="border_length", edge_perimeter_col="length")
    new_graph = MultiLevelGraph(new_base_graph, ["block_name"]);
    push!(TH.graphs_by_level, new_graph)
end

function delete_level_at_top!(TH)
    pop!(TH.graphs_by_level)
    pop!(TH.children_by_level)
    pop!(TH.parents_by_level) 
    pop!(TH.onbrs_by_level)
end

function manually_undo_merge_at_top_level(TH, node)


end


function TemperingHierarchy(Gbase, filename_prefix, num_levels::Int)
    graphs_by_level = [Gbase]
    children_by_level = [Dict{GlobalID, Vector{GlobalID}}()]
    parents_by_level = []
    for i=2:num_levels 
        #node_data = Set(["block_name", "pop2020cen", "area", "border_length", "G16_PR_D", "G16_PR_R", "id"])
        #new_base_graph = BaseGraph(filename_prefix*string(i)*".json", "pop2020cen", inc_node_data=node_data, area_col="area",node_border_col="border_length", edge_perimeter_col="length")
        node_data = Set(["block_name", Gbase.base_graph.pop_col, "area", "border_length", "id"])
        new_base_graph = BaseGraph(filename_prefix*string(i)*".json", Gbase.base_graph.pop_col, inc_node_data=node_data, area_col="area",node_border_col="border_length", edge_perimeter_col="length")
        new_graph = MultiLevelGraph(new_base_graph, ["block_name"]);
        push!(graphs_by_level, new_graph)
    end

    for i=2:num_levels 
        d = JSON.parsefile(filename_prefix*string(i)*"_children.json")
        d2 = Dict{GlobalID, Vector{GlobalID}}([parse(Int, k) => v for (k,v) in d])
        push!(children_by_level, d2)

        d = JSON.parsefile(filename_prefix*string(i-1)*"_parents.json")
        d2 = Dict{GlobalID, GlobalID}([parse(Int, k) => v for (k,v) in d])
        push!(parents_by_level, d2)
    end
    return TemperingHierarchy(graphs_by_level, children_by_level, parents_by_level, [], filename_prefix)
end

function highest_ancestor(TH, node_ID) 
    return ancestor_at_level(TH, node_ID, length(TH.graphs_by_level))    
end

function ancestor_at_level(TH, node_ID, level, cur_level=1)
    @assert node_ID in finest_level_nodes(TH.graphs_by_level[cur_level]; is_flat=true)
    if cur_level == level
        return node_ID 
    end 
    return ancestor_at_level(TH, TH.parents_by_level[cur_level][node_ID], level, cur_level+1)
end


function onbrs_if_matching_made(Gnc, matching, onbrs) 
    all_nodes = finest_level_nodes(Gnc) 
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

    node_ID = maximum(finest_level_nodes(Gnc)) + 1000
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

function remove_artpt_forming_edges(matching, G, og_onbrs)
    matching_to_use = []
    c_onbrs, c_pair_parent = onbrs_if_matching_made(G, matching, og_onbrs) 
    for edge in matching 
        if count(x -> x == edge.id1, og_onbrs[edge.id2]) > 1 || 
            count(x -> x == edge.id2, og_onbrs[edge.id1]) > 1 
            @assert count(x -> x==edge.id1, og_onbrs[edge.id2]) > 1
            @assert count(x -> x==edge.id2, og_onbrs[edge.id1]) > 1
            println("encloses a space")
            continue
        elseif count(x -> x == -1, c_onbrs[c_pair_parent[edge.id1]]) > 1 
            println("touches edge twice") 
            continue  
        else 
            push!(matching_to_use, edge)
        end
    end
    return matching_to_use
end

# guarantees that nodes that are NOT paired (i.e., finest level 
# nodes of Gnc that are not elements of matching) have the SAME node ID in 
# the parent graph as in Gnc. 
function pair_up_nodes(Gnc, matching, output_file, new_level_num, onbrs, prefix="pct")
    blocks_dict = Dict() 
    blocks_dict["nodes"] = []
    these_nodes = [] 

    matched_nodes = Set()
    for edge in matching 
        push!(matched_nodes, edge.id1) 
        push!(matched_nodes, edge.id2)
    end

    pair_parent = Dict{GlobalID, GlobalID}()
    pair_children = Dict{GlobalID, Vector{GlobalID}}()
    new_onbrs = Dict{GlobalID, Vector{Vector{GlobalID}}}()

    all_nodes = finest_level_nodes(Gnc) 
    for n1 in all_nodes 
        if n1 in matched_nodes continue end 
        node_ID = n1 
        pair_parent[n1] = node_ID
        pair_children[node_ID] = [n1]
        push!(these_nodes, node_ID)
        this_block_info = Dict() 
        this_block_info["id"] = node_ID
        this_block_info["block_name"] = prefix*"single_"*string(node_ID)
        this_block_info[Gnc.base_graph.pop_col] = get_node_population(Gnc, n1)
        this_block_info["area"] = node_attributes(Gnc, n1)["area"] 
        this_block_info["border_length"] = node_attributes(Gnc, n1)["border_length"]
        #this_block_info["G16_PR_D"] = node_attributes(Gnc, n1)["G16_PR_D"] 
        #this_block_info["G16_PR_R"] = node_attributes(Gnc, n1)["G16_PR_R"]
        push!(blocks_dict["nodes"], this_block_info)
        new_onbrs[n1] = [onbrs[n1]]
    end

    node_ID = maximum(finest_level_nodes(Gnc)) + 1000
    for edge in matching
        node_ID += 1 
        n1, n2 = edge.id1, edge.id2
        pair_parent[n1] = node_ID
        pair_parent[n2] = node_ID
        pair_children[node_ID] = [n1, n2]
        push!(these_nodes, node_ID)
        this_block_info = Dict() 
        this_block_info["id"] = node_ID
        this_block_info["block_name"] = prefix*"pair_"*string(node_ID)
        this_block_info[Gnc.base_graph.pop_col] = get_node_population(Gnc, n1) + get_node_population(Gnc, n2)
        this_block_info["area"] = node_attributes(Gnc, n1)["area"] + node_attributes(Gnc, n2)["area"]
        this_block_info["border_length"] = node_attributes(Gnc, n1)["border_length"] + node_attributes(Gnc, n2)["border_length"]
        #this_block_info["G16_PR_D"] = node_attributes(Gnc, n1)["G16_PR_D"] + node_attributes(Gnc, n2)["G16_PR_D"]
        #this_block_info["G16_PR_R"] = node_attributes(Gnc, n1)["G16_PR_R"] + node_attributes(Gnc, n2)["G16_PR_R"]
        push!(blocks_dict["nodes"], this_block_info)
        new_onbrs[node_ID] = combine_onbrs(n1, n2, onbrs)
    end

    # finish dealing with onbrs 
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

    blocks_dict["adjacency"] = []

    edge_lens = Dict{UndirectedEdge, Real}()
    for edge in finest_level_edges(Gnc) 
        p1, p2 = pair_parent[edge.id1], pair_parent[edge.id2]
        if p1 == p2 continue end 
        pedge = UndirectedEdge(p1, p2)
        if haskey(edge_lens, pedge) 
            edge_lens[pedge] += get_edge_attributes(Gnc, edge).length 
        else 
            edge_lens[pedge] = get_edge_attributes(Gnc, edge).length 
        end
    end

    for node in these_nodes 
        adj_pairs = Set([pair_parent[x] for y in pair_children[node] for (x, _) in get_finest_level_neighbors(Gnc, y)]) 
        this_adj = []
        for nbr_ID in adj_pairs
            if nbr_ID == node continue end 
            push!(this_adj, Dict(["id" => nbr_ID, "length" => edge_lens[UndirectedEdge(node, nbr_ID)]]))
        end
        push!(blocks_dict["adjacency"], [node, this_adj])
    end

    blocks_dict["directed"] = false
    blocks_dict["multigraph"] = false
    blocks_dict["graph"] = []

    open(output_file*string(new_level_num)*".json","w") do f
        JSON.print(f, blocks_dict,4)
    end
    open(output_file*string(new_level_num-1)*"_parents.json","w") do f
        JSON.print(f, pair_parent, 4)
    end
    open(output_file*string(new_level_num)*"_children.json","w") do f
        JSON.print(f, pair_children, 4)
    end

    return pair_parent, pair_children, coarsened_onbrs
end


function write_TH_partition_to_file(TH, partition, lvl, output_file)
    @assert partition.graph == TH.graphs_by_level[lvl]
    G1 = TH.graphs_by_level[1]
    out = []
    for node_ID in finest_level_nodes(G1; is_flat=true)
        dist = partition.node_to_district[ancestor_at_level(TH, node_ID, lvl)]
        node_name = node_name_exclude_quotient_nodes(G1, node_ID)
        s = join(node_name, ",")
        push!(out, "\""*s*"\",\""*string(dist)*"\"")
    end
    s = "node_name,district\n"*join(out, "\n")
    open(output_file, "w") do file 
        write(file, s)
    end 
end

function write_cpartition_to_file(cpartition, pair_children, Gncpct, output_file) 
    out = []
    for (k,v) in cpartition.node_to_district
        if k == cpartition.graph.root.global_ID continue end  
        for c in pair_children[k]
            push!(out, "\""*s*"\",\""*string(v)*"\"")
        end
    end
    
end