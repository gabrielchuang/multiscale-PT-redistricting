# returns (d1, d2, d3), p, where p=1/num such subgraphs 
# unless just_count=true, in which case just returns p 
function uniformly_pick_size_3_connected_induced_subgraph(partition, rng; just_count = false)
    possibilities = [] 
    num_possibilities = 0 
    for i=1:partition.num_dists-2 
        for j=i+1:partition.num_dists-1 
            for k=j+1:partition.num_dists 
                if count([partition.adjacency[i,j], 
                    partition.adjacency[j,k], partition.adjacency[i,k]]) >= 2
                    if just_count 
                        num_possibilities += 1 
                    else 
                        push!(possibilities, (i,j,k))
                    end 
                end 
            end
        end
    end
    if just_count
        return 1.0 / num_possibilities
    else 
        return uniformly_choose(possibilities, rng)
    end
end 

function multistep_forest_recom2!(
    partition::MultiLevelPartition, 
    measure::Measure,
    constraints::Dict,
    rng::AbstractRNG, 
    i::Int;
    step_size=10
) 
    steps_taken = 0 
    max_attempts = 100
    num_zero_probs = 0 

    districts_to_choose_from, prob_region_cur = uniformly_pick_size_3_connected_induced_subgraph(partition, rng)
    districts_to_choose_from = collect(districts_to_choose_from)

    # if we end up rejecting, then we'll have to undo all the moves. 
    snapshot_og = take_snapshot(partition, districts_to_choose_from) 
    function revert_multimove() 
        restore_from_snapshot!(partition, districts_to_choose_from, snapshot_og...)
    end

    # initialized to 1; will update as we make steps 
    overall_transition_prob = 1
    for attempt = 1:max_attempts
        # pick a random district 
        
        measure_g0 = Measure(0)
        p, update = forest_recom2_2!(partition, measure_g0, constraints, rng, i; 
            districts_to_choose_from=districts_to_choose_from) 

        if update === nothing || p == 0 
            num_zero_probs += 1 
        else 
            overall_transition_prob *= p
            update_partition!(partition, update) 
            steps_taken += 1
        end 

        if steps_taken == step_size 
            push!(partition.extensions[MULTISTEP_ATTEMPTS], attempt)

            snapshot_new = take_snapshot(partition, districts_to_choose_from) 
            function accept_multimove() 
                restore_from_snapshot!(partition, districts_to_choose_from, snapshot_new...)
            end

            log_l_edge_prob_proposed = log_linking_edge_prob_tree_space(partition)
            log_tree_count_proposed = log_hierarchical_forest_count(partition)
            prob_region_proposed = uniformly_pick_size_3_connected_induced_subgraph(partition, rng, just_count=true)

            delta_energy = 0 
            for ii = 1:length(measure.weights)
                weight = measure.weights[ii]
                if weight != 0
                    energy = measure.scores[ii]
                    delta_energy -= weight*energy(partition, districts_to_choose_from)
                end
            end

            revert_multimove()
        
            for ii = 1:length(measure.weights)
                weight = measure.weights[ii]
                if weight != 0
                    energy = measure.scores[ii]
                    delta_energy += weight*energy(partition, districts_to_choose_from)
                end
            end

            log_l_edge_prob_cur = log_linking_edge_prob_tree_space(partition)
            log_tree_count_cur = log_hierarchical_forest_count(partition)

            log_linking_edge_ratio = log_l_edge_prob_proposed - log_l_edge_prob_cur
            log_tree_count_ratio = log_tree_count_proposed - log_tree_count_cur 
            prob_region_ratio = prob_region_cur / prob_region_proposed

            pi_ratio = exp(measure.gamma*(log_linking_edge_ratio + log_tree_count_ratio))

            push!(partition.extensions[MS_QRATIO], overall_transition_prob)
            push!(partition.extensions[MS_REGIONRATIO], prob_region_ratio)
            push!(partition.extensions[MS_PIRATIO], pi_ratio)
            push!(partition.extensions[MS_EXP_DELTAENERGY], exp(delta_energy))

            overall_transition_prob *= prob_region_ratio 
            overall_transition_prob *= pi_ratio 
            overall_transition_prob *= exp(delta_energy)
        
            return overall_transition_prob, (revert_multimove, accept_multimove, num_zero_probs) 
        end
    end

    # we failed max_attempts times; revert it and return prob 0 
    revert_multimove() 
    return 0, (nothing, nothing, num_zero_probs)
end


function build_multistep_forest_recom2(
    constraints::Dict; 
    step_size::Int=10
)
    f(p, m, r, i) = multistep_forest_recom2!(p, m, constraints, r, i, step_size=step_size)
    return f
end
