
function get_vap(G::MultiLevelGraph)::Real 
    node_attributes(G, G.root.global_ID)[G.base_graph.vap_col]
end

function get_bvap(G::MultiLevelGraph)::Real 
    node_attributes(G, G.root.global_ID)[G.base_graph.bvap_col]
end

function get_bpop(G::MultiLevelGraph)::Real 
    node_attributes(G, G.root.global_ID)[G.base_graph.bpop_col]
end

function get_vra_scores(
    partition::MultiLevelPartition,
    districts = 1:partition.num_dists
)
    pop = partition.graph.base_graph.total_pop
    bpop = get_bpop(partition.graph)
    num_dists = partition.num_dists
    target_mino_dists = Int(floor(num_dists*bpop/pop))
    # @show num_dists, bpop/pop, num_dists*bpop/pop, target_mino_dists

    vapfracs = Vector{Real}(undef, num_dists)

    for di = 1:partition.num_dists
        vap = get_vap(partition.subgraphs[di])
        bvap = get_bvap(partition.subgraphs[di])
        vapfracs[di] = bvap/vap
    end
    sort!(vapfracs, rev=true)
    return vapfracs[1:target_mino_dists]
end

function get_vra_score(
    partition::MultiLevelPartition,
    districts = 1:partition.num_dists
)
    score = 0
    vapfracs = get_vra_scores(partition)

    for ii = 1:length(vapfracs)
        if vapfracs[ii] < 0.5
            score += sqrt(0.5-vapfracs[ii])
        end
    end
    return score
end