@enum EnergyType POPULATION COUNTY_SPLITS COMPACTNESS COUNTY_SPLITS_MULTIPLICITY COUNTY_SPLITS_DELTA CUT_EDGES COMPACTNESS_EXP PATH_COMPACTNESS FLIPS EIGHTD_SPLITS

# the measure is e^-J(x)
# so, J(x) is something like J_pop(x) + J_split(x) + J_compact(x) 
# J_pop(x) = w_pop (sum_d (pop(d) - ideal)^2) * 1(|pop_d - ideal| > max)
# J_split(x) = w_split * (# splits over max) 
# J_compact(x) = w_compact * sum_d iso_score(d)

# We would like to be able to partially compute the energy (e.g., compute 
# the compactness of each district, then hand it over to the caller to potentially 
# make changes (e.g., if SNF, we can change the compactness of just those two dists) 
# and then finish computing it (e.g., by summing up the iso scores)). 
# Allows for overall constant computation of energies for individual SNF proposals 
# (important for tempered proposals) 
struct Energy 
    full_compute::Function # partition -> Real
    partial_compute::Function # partition -> T 
    finish_compute::Function # T -> Real 
    energy_type::EnergyType
end

# J_pop(x) = w_pop (sum_d (pop(d) - ideal) * 1(|pop_d - ideal| > max) )
# partial is < pop_i for i in dists> 
# finish_compute is (sum_d (pop(d) - ideal) * 1(|pop_d - ideal| > max) )
function populationEnergy(ideal_pop::Real, allowable_percent_deviation::Float64)::Energy
    function partial_compute(
        partition::MultiLevelPartition; 
        dists=1:partition.num_dists
    )::Dict{Int, Int} 
        res = Dict{Int, Int}()
        for dist in dists 
            res[dist] = partition.dist_populations[dist]
        end
        return res 
    end

    function finish_compute(partial::Dict{Int, Int})
        res = 0 
        for (k, dist_pop) in partial 
            if abs(dist_pop - ideal_pop) > allowable_percent_deviation * ideal_pop 
                res += ((abs(dist_pop - ideal_pop) / ideal_pop))
            end
        end
        return res
    end

    function full_compute(partition::MultiLevelPartition; dists=1:partition.num_dists)
        return finish_compute(partial_compute(partition; dists=dists)) 
    end
    
    return Energy(full_compute, partial_compute, finish_compute, POPULATION)
end

# J_compact(x) = w_compact * sum_d iso_score(d)
# partial is <area(d1) ... area(dn)>, <perim(d1), ... perim(dn)>
# finish_compute is sum_di perim(di)^2 / area(di) 
function compactnessEnergy(reweight)::Energy 
    function partial_compute(
        partition::MultiLevelPartition;
        dists=1:partition.num_dists
    )::Tuple{Dict{Int, Float64}, Dict{Int, Float64}}
        areas, perims = Dict{Int, Float64}(), Dict{Int, Float64}() 
        for di in dists 
            areas[di] = partition.dist_areas[di] 
            perims[di] = partition.dist_perimeters[di]
            #isoperimetric_score = pr^2/area
        end
        return (areas, perims)
    end

    function finish_compute(partial::Tuple{Dict{Int, Float64}, Dict{Int, Float64}})
        areas, perims = partial 
        res = 0
        isos = sort([perims[k] for k=1:length(perims)] .^2 ./ [areas[k] for k=1:length(areas)])
        
		if reweight !== nothing
			return sum(isos) - 0.7 * isos[end] - 0.5 * isos[end-1]
		end 
		return sum(isos) 

        #=for i in keys(areas)
            res += (perims[i] ^ 2 / areas[i]) ^ norm 
        end
        return res ^ (1/norm)=#
    end

    function full_compute(partition::MultiLevelPartition; dists=1:partition.num_dists)
        return finish_compute(partial_compute(partition; dists=dists)) 
    end
    
    return Energy(full_compute, partial_compute, finish_compute, COMPACTNESS)
end 

function articulationEnergy(threshold)::Energy 
    function partial_compute(partition, dists) 
        @assert false "can't use articulation energy to temper" 
    end 
    function finish_compute(partial) 
        @assert false "can't use articulation energy to temper"
    end 
    function full_compute(partition::MultiLevelPartition; dists=1:partition.num_dists) 
        return maximum([0, sum([length(x) for (_, x) in partition.articulation_points])-threshold])
    end
    return Energy(full_compute, partial_compute, finish_compute, COMPACTNESS)
end


function compactnessEnergyExperimental()::Energy 
    function partial_compute(
        partition::MultiLevelPartition;
        dists=1:partition.num_dists
    )::Tuple{Dict{Int, Float64}, Dict{Int, Float64}, Float64}
        areas, perims = Dict{Int, Float64}(), Dict{Int, Float64}() 
        og_compactness = 0
        for di in dists 
            areas[di] = partition.dist_areas[di] 
            perims[di] = partition.dist_perimeters[di]
            #isoperimetric_score = pr^2/area
            og_compactness += perims[di]^2 / areas[di]
        end
        return (areas, perims, og_compactness)
    end

    function finish_compute(partial::Tuple{Dict{Int, Float64}, Dict{Int, Float64}, Float64})
        areas, perims, og_compactness = partial 
        res = 0
        for i in keys(areas)
            res += perims[i]^2 / areas[i]
        end
        if res < og_compactness # more compact 
            return res / 2.0
        else 
            return res 
        end
    end

    function full_compute(partition::MultiLevelPartition; dists=1:partition.num_dists)
        return finish_compute(partial_compute(partition; dists=dists)) 
    end
    
    return Energy(full_compute, partial_compute, finish_compute, COMPACTNESS_EXP)
end 

function splitNodesEnergy(Gnc, Gnc8d, max_splits)::Energy
    function partial_compute(partition::MultiLevelPartition; dists=1:partition.num_dists)
        return length(split_nodes(Gnc, Gnc8d, partition)), Gnc
    end

    function finish_compute(partial)
        return maximum([0, partial[1] - max_splits])
    end

    function full_compute(partition::MultiLevelPartition; dists=1:partition.num_dists)
        return finish_compute(partial_compute(partition)) 
    end
    
    return Energy(full_compute, partial_compute, finish_compute, EIGHTD_SPLITS)
end 




# J_split(x) = w_split * (# splits over max) 
# the intermediate (partially-computed) form is the number of split counties. 
# the output energy is # split counties - max_splits 
function count_split_counties(partition::MultiLevelPartition) 
    split_counties = 0
    for (node, dist) in partition.node_to_district 
        name = node_name(partition.graph, node)
        if dist == -1 && name !== nothing && length(name) == 2
            split_counties += 1 
        end
    end
    return split_counties 
end


function countySplitsEnergy(max_splits::Int)::Energy 
    function partial_compute(partition::MultiLevelPartition)::Int
        return count_split_counties(partition)
    end

    function finish_compute(partial::Int)
        return maximum([0, partial - max_splits])
    end

    function full_compute(partition::MultiLevelPartition)
        return finish_compute(partial_compute(partition)) 
    end
    
    return Energy(full_compute, partial_compute, finish_compute, COUNTY_SPLITS)
end 

function districts_in_county(partition::MultiLevelPartition, cty::GlobalID) 
    n2d = partition.node_to_district
    @assert haskey(n2d, cty) 
    if n2d[cty] == -1 
        return Set([n2d[x] for x in get_children_IDs(partition.graph, cty)])
    else 
        return Set([n2d[cty]])
    end
end


# J_split(x) = w_split * (# splits over max) 
# the intermediate (partially-computed) form is the number of split counties. 
# the output energy is # split counties - max_splits 
function count_split_counties_multiplicity(partition)
    split_counties = 0
    for (node, dist) in partition.node_to_district 
        if dist == -1 && partition.graph.distance_from_root[node] == 1
            split_counties += (length(districts_in_county(partition, node))-1)^2
        end
    end
    return split_counties 
end

function countySplitsMultiplicityEnergy(max_splits::Int)::Energy 
    function partial_compute(partition::MultiLevelPartition; dists=1:partition.num_dists)::Int
        count_split_counties_multiplicity(partition)
    end

    function finish_compute(partial::Int)
        return maximum([0, partial - max_splits])
    end

    function full_compute(partition::MultiLevelPartition; dists=1:partition.num_dists)
        return finish_compute(partial_compute(partition; dists=dists)) 
    end
    
    return Energy(full_compute, partial_compute, finish_compute, COUNTY_SPLITS_MULTIPLICITY)
end 

function deltaSplitsEnergy()::Energy 
    function partial_compute(partition; dists=[]) 
        return 0 
    end
    function finish_compute(partial::Int) 
        return partial 
    end 
    function full_compute(partition::MultiLevelPartition; dists=[])
        return finish_compute(partial_compute(partition; dists=dists)) 
    end
    
    return Energy(full_compute, partial_compute, finish_compute, COUNTY_SPLITS_DELTA)
end 
    

function cutEdgeEnergy()::Energy 
    function partial_compute(partition; dists=[]) 
        return count_finest_level_cut_edges(partition)
    end 
    function finish_compute(partial::Int) 
        return partial 
    end 
    function full_compute(partition::MultiLevelPartition; dists=[])
        return finish_compute(partial_compute(partition; dists=dists)) 
    end
    
    return Energy(full_compute, partial_compute, finish_compute, CUT_EDGES)
end 


function leaf_nodes_of(partition, dist) 
    res = Set{GlobalID}()
    for (k,v) in partition.district_to_nodes[dist]
        if partition.graph.distance_from_root[k] == 1 
            if v === nothing 
                union!(res, get_children_IDs(partition.graph, k))
            else 
                union!(res, v) 
            end
        end
    end
    return res 
end 

function pathCompactness()::Energy 
    function partial_compute(partition; dists=[]) 
        res = Dict{Int, Dict{GlobalID, Float64}}()
        for dist in dists 
            dist_res = Dict{GlobalID, Float64}()

            dist_pcts = collect(leaf_nodes_of(partition, dist))
            g, vmap = induced_subgraph(partition.graph.base_graph.simple_graph, dist_pcts)
            trav = floyd_warshall_shortest_paths(g) 
            pts = enumerate_paths(trav)
            for i in eachindex(pts) 
                avg_path_length = sum(length(pts[i])) / length(pts[i])
                dist_res[vmap[i]] = avg_path_length 
            end 
            avg_avg = sum([v for (_,v) in dist_res]) / length(dist_res) 
            for (k,v) in dist_res 
                dist_res[k] = avg_avg / v
            end 
            res[dist] = dist_res 
        end
        return res 
    end
    function finish_compute(partial)
        return 0
    end 
    function full_compute(partition::MultiLevelPartition; dists=[])
        return finish_compute(partial_compute(partition; dists=dists)) 
    end
    
    return Energy(full_compute, partial_compute, finish_compute, PATH_COMPACTNESS)
end

# weight the potential flips by the shortest path to the *next* district in the cycle. 
# for example, if we have d1 -> d2 -> d3 (-> d1), and a pct in d2 being eaten by d1 is weighted 
# by its distance to d3. 
function cyclicPathCompactness()::Energy 
    function partial_compute(partition; dists=[]) 
        res = Dict{Int, Dict{GlobalID, Float64}}()

        for i in eachindex(dists)
            iminus1 = i == 1 ? length(dists) : i-1 

            dist = dists[i]

            dist_res = Dict{GlobalID, Float64}()            
            dist_pcts = collect(leaf_nodes_of(partition, dist))
            g, vmap = induced_subgraph(partition.graph.base_graph.simple_graph, dist_pcts)

            artificial_node = nv(g) + 1 # we'll use this node to represent the next district (e.g. d3)
            add_vertex!(g)
            push!(vmap, -1)
            @assert artificial_node in vertices(g) 

            for node_idx in eachindex(vmap) 
                if vmap[node_idx] in partition.adjacent_fine_neighbors[dists[iminus1], dist]
                    # it's on the border; add an edge 
                    add_edge!(g, node_idx, artificial_node)
                end
            end

            trav = dijkstra_shortest_paths(g, artificial_node).dists
            avg_length = sum(trav) / length(trav) 
            stdev_length = std(trav)

            z_score = i -> trav[i] - avg_length / stdev_length
            pos_z_score = i -> z_score(i) > 0 ? z_score(i) : 0 

            dist_res = Dict{GlobalID, Float64}([vmap[i] => 10-pos_z_score(i)-1 for i in eachindex(trav)])
            res[dist] = dist_res 
        end
        return res 
    end
    function finish_compute(partial)
        return 0
    end 
    function full_compute(partition::MultiLevelPartition; dists=[])
        return 0#finish_compute(partial_compute(partition; dists=dists)) 
    end
    
    return Energy(full_compute, partial_compute, finish_compute, PATH_COMPACTNESS)
end

function num_valid_flips_energy()::Energy 
    function partial_compute(partition; dists=[]) 
        return partition, dists        
    end
    function finish_compute(partial)
        partition, dists = partial
        res = 0
        ap = partition.articulation_points

        for i in eachindex(dists)
            iminus1 = i == 1 ? length(dists) : i-1 

            npcs1 = Set(
                filter(x -> !(x in ap[dists[iminus1]]), 
                collect(partition.adjacent_fine_neighbors[dists[i], dists[iminus1]]))
            )
            npcs2 = Set(
                filter(x -> !(x in ap[dists[i]]), 
                collect(partition.adjacent_fine_neighbors[dists[iminus1], dists[i]]))
            )
            if length(npcs1) >= length(npcs2)
                res += 0.5*(length(npcs1) - length(npcs2))
            else 
                res += (length(npcs2) - length(npcs1))^2
            end 
        end
        return res 
    end 
    function full_compute(partition::MultiLevelPartition; dists=[])
        return finish_compute(partial_compute(partition; dists=dists)) 
    end

    return Energy(full_compute, partial_compute, finish_compute, FLIPS)
end


function spanningForestEnergy(gamma)::Energy 
    function partial_compute(partition) 
        @assert false "can't temper with spanning forest energy"
    end 
    function finish_compute(partial) 
        @assert false "can't temper with spanning forest energy"
    end 
    function full_compute(partition::MultiLevelPartition; dists=[])
        if dists != [] && length(dists) != partition.num_dists
            @assert false "hmm"
        end
        return (1-gamma) * get_log_spanning_forests(partition) 
    end
    
    return Energy(full_compute, partial_compute, finish_compute, CUT_EDGES)
end 
