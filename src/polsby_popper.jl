
function get_area(G::MultiLevelGraph)::Real 
    node_attributes(G, G.root.global_ID)[G.base_graph.area_col]
end

function traverse_perimeter(partition::MultiLevelPartition, nbr::MultiLevelNeighbor, d; 
    pshow=false
) 
    @assert nbr.global_ID != 0
    dist = look_up_district(partition, nbr.global_ID) 
    if dist == -1 
        return sum([traverse_perimeter(partition, x, d) for x in values(nbr.finer_neighbors)])
    elseif dist == d 
        return 0 
    else 
        if pshow 
            println(node_name(partition.graph, nbr.global_ID), " ", nbr.edge_attributes.length)
        end 
        return nbr.edge_attributes.length 
    end 
end 

function get_perimeter(partition::MultiLevelPartition, d::Int)::Real
    d2n = partition.district_to_nodes[d]
    G = partition.graph 

    border_length = 0
    for (k,v) in d2n 
        if v === nothing 
            for nbr in fine_neighbors(partition.graph, k) 
                border_length += traverse_perimeter(partition, nbr, d)
            end
        end
    end

    external_border_length = node_attributes(partition.subgraphs[d], G.root.global_ID)[G.base_graph.node_border_col]
    return border_length + external_border_length
end

function get_isoperimetric_score(
    partition::MultiLevelPartition,
    districts = 1:partition.num_dists
)
    score = 0
    for di in districts
        area = get_area(partition.subgraphs[di])
        pr = get_perimeter(partition, di)
        isoperimetric_score = pr^2/area
        score += isoperimetric_score
    end
    return score
end

#=
function get_perimeter(
    graph::MultiLevelGraph, 
    flat_district_leaves 
)
    border_length = 0 
    for node in flat_district_leaves
        for nbr in fine_neighbors(graph, node) 
            if !(nbr.global_ID in flat_district_leaves) 
                border_length += nbr.edge_attributes.length
            end 
        end
    end 
    return border_length 
end=#