abstract type AbstractConstraint end

# @enum CONSTRAINT pop pack_node contiguous_traversal max_coarse_splt max_shared_coarse excess_dists_in_coarse max_hamming_distance max_shared_nodes

struct PopulationConstraint <: AbstractConstraint
    min_pop::Real
    max_pop::Real
end

struct ConstrainDiscontinuousTraversals <: AbstractConstraint
    max_line_segments::Int
end

struct MaxSplitNodes <: AbstractConstraint 
    max_split_nodes::Function # from MultiLevelGraph * GlobalID to Int. a given node may not have more than max(x) of its 
    # children split between districts. 
end 

struct MaxSplitNodesByLevel <: AbstractConstraint 
    max_split_nodes::Int 
end 
