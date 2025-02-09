mutable struct Measure
    gamma::Float64
    weights::Vector{Float64}
    scores::Vector{Energy}
    descriptions::Vector{String}
end

# partially computes each energy: generates a length-k vector if there are k energies
function partially_compute_energy(
    partition::MultiLevelPartition, 
    measure::Measure; 
    dists=1:partition.num_dists
)::Vector{Any} 
    partially_computed_energies = Vector{Any}() 
    for score in measure.scores
        push!(partially_computed_energies, score.partial_compute(partition; dists=dists))
    end
    return partially_computed_energies
end

# given the partially computed energies, finishes the computation 
function finish_compute_energy(
    measure::Measure, 
    partially_computed_energies::Vector{Any} 
)::Float64 
    score = 0.0
    for i=1:length(partially_computed_energies)
        score += measure.weights[i] * measure.scores[i].finish_compute(partially_computed_energies[i])
    end 
    return score 
end

function fully_compute_energy(partition::MultiLevelPartition, measure::Measure, dists=1:partition.num_dists) 
    res = 0
    for i in eachindex(measure.scores)
        res += measure.weights[i] * measure.scores[i].full_compute(partition; dists=dists)
    end
    return res
end 

""""""
function Measure(gamma::Float64)
    scores = Vector{Energy}(undef, 0)
    weights = Vector{Float64}(undef, 0)
    descriptions = Vector{String}(undef, 0)
    return Measure(gamma, weights, scores, descriptions)
end

function push_energy!(measure::Measure, energy::Energy, weight::Real, desc::String) 
    push!(measure.scores, energy) 
    push!(measure.weights, weight) 
    push!(measure.descriptions, desc) 
end

function get_delta_energy(
    partition::MultiLevelPartition, 
    measure::Measure,
    update::Tuple
)
    score = 0

    changed_districts, new_dist_masks, new_dist_pops, edge = update
    for ii = 1:length(measure.weights)
        weight = measure.weights[ii]
        if weight == 0
            continue
        end
        energy = measure.scores[ii]
        score += weight*energy.full_compute(partition; dists=changed_districts)
    end

    snapshot = take_snapshot_and_update!(partition, changed_districts, new_dist_masks, new_dist_pops) 

    for ii = 1:length(measure.weights)
        weight = measure.weights[ii]
        if weight == 0
            continue
        end
        energy = measure.scores[ii]
        score -= weight*energy.full_compute(partition; dists=changed_districts)
    end
   
    restore_from_snapshot!(partition, snapshot...)

    return exp(score)
end

function log_measure(
    partition::MultiLevelPartition,
    measure::Measure;
    gamma::Real=measure.gamma
)
    log_linking_edges = get_log_linking_edges(partition)
    log_forests = get_log_spanning_forests(partition)
    log_p = (1-gamma)*(log_forests + log_linking_edges) -
        fully_compute_energy(partition, measure)

    return log_p
end

function get_log_spanning_forests(partition::MultiLevelPartition)
    # requires a flat graph
    lnsp = 0
    for dist_subgraph in partition.subgraphs
        s_graph = child_graph(dist_subgraph, dist_subgraph.root.global_ID)[1]
        lnsp += log_nspanning(s_graph)
    end
    return lnsp 
end
