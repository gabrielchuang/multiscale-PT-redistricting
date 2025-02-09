module MultiScaleMapSampler
using JSON
using SimpleWeightedGraphs
using Graphs
using RandomNumbers
using LinearAlgebra
using Hungarian
using Printf
using Combinatorics
using DataStructures
using Statistics
using InvertedIndices 

export AbstractGraph,
    BaseGraph,

    # graph functions 
    MultiLevelGraph, 
    MultiLevelPartition,
    fine_neighbors, 
    take_subgraph,
    get_children_IDs, 
    get_edge_attributes,
    child_graph,
    node_attributes, 
    MultiLevelNeighbor, 
    MultiLevelNode, 
    GlobalID,
    EdgeAttributes,
    log_hierarchical_tree_count,
    precompute_log_tree_counts!,
    log_hierarchical_forest_count, 
    node_name, is_finest_level_node,

    # forest recom mcmc 
    run_metropolis_hastings!, 
    build_forest_recom2, 
    build_multistep_forest_recom2_2, 

    # partition utilities 
    write_partition_to_file,
    look_up_district,
    write_graph_to_file,
    write_hierarchy_to_file, 
    print_dict, 
    any_pct_in,
    dists_adj_to_node,
    partition_from_d2n,
    output,

    # graph utilities 
    is_intact,
    get_edge_weight,
    is_in_graph,
    get_children_IDs, 
    compute_all_fine_neighbors,
    preord_foreach, 
    is_descendant, 
    ancestor_at_depth, 
    get_node_population, 
    get_finest_level_neighbors,
    get_parent_ID,
    write_graph_nodes_to_file,
    node_name_exclude_quotient_nodes,
    UndirectedEdge, 
    is_in_graph,
    finest_level_edges,
    finest_level_nodes,
    finest_level_edges_adjlist,
    get_node_area,
    filter_edges!,
    get_node_perimeter,
    graph_subset_by_cty,

    # single node flip and multistep single node flip
    multistep_SNF, 
    construct_county_split_graph, 
    single_node_flip2,
    run_metropolis_hastings_SNF!,
    accept_SNF!,

    REJECTION_COUNT, ACCEPTANCE_COUNT, TRANSITION_RATIOS,
    TRANSITION_RATIOS2,
    
    # masks
    mask_intersection,
    mask_union,  
    clean_mask!, 

    # dynamic hierarchies 
    expand_hierarchy!, 

    # constraints
    initialize_constraints,
    add_constraint!,
    AbstractConstraint,
    PopulationConstraint,
    ConstrainDiscontinuousTraversals,
    PackNodeConstraint,
    MaxCoarseNodeSplits,
    MaxSharedCoarseNodes,
    MaxSharedNodes,
    AllowedExcessDistsInCoarseNodes,
    MaxHammingDistance,
    MaxSplitNodes, 
    MaxSplitNodesByLevel,

    # Writer
    Writer, 
    close_writer,
    write_districts, 
    SimpleWriter, 
    push_writer!,

    # energies/measure
    Measure,
    push_measure!,
    push_energy!,
    populationEnergy, 
    compactnessEnergy, 
    countySplitsEnergy, 
    fourthPowerCompactnessEnergy,
    countySplitsMultiplicityEnergy,
    deltaSplitsEnergy,
    cutEdgeEnergy,
    compactnessEnergyExperimental,
    pathCompactness,
    cyclicPathCompactness,
    fully_compute_energy,
    count_finest_level_cut_edges,
    get_log_spanning_forests,
    get_log_spanning_trees,
    spanningForestEnergy,

    get_isoperimetric_score,
    get_area, 
    get_perimeter,

    # plane graphs / Quad Edge / perfect matchings 
    QuadEdgeMesh, 
    count_perfect_matchings,
    log_count_perfect_matchings,
    construct_pfaffian_orientation, 
    sample_perfect_matching, 
    fast_sample_pm,
    UndirectedEdge,
    even_ccs,
    symmetrize_0s_old!,
    symmetrize_0s_faster!,
    add_level_with_pop_threshold!,
    filter_edges!,
    
    # hierarchical partition 
    combine_onbrs, 
    multisplice_at, 
    splice_at, 
    remove_adj_duplicates_cyclic!,
    add_level_at_top!,
    remove_level_at_top!,
    node_pop,
    add_level_with_pop_target!,
    remove_multiedges!,
    HierarchicalPartition,
    topmost_ancestor,
    add_level_with_edge_reweight!,
    isoperimetric_score,
    add_descendants_at_level!,
    spanning_forest_count,
    pm_count,
    iso_matching,
    g_ancestor,
    matching_maximize_below_thres,
    iso_matching2,
    iso_matching3,

    # parallel tempering 
    parallel_tempering!, 
    Chain, 
    run_chain!, 
    try_swap_replicas!,
    get_vra_score,

    # purely for testing 
    enumerate_triangles, 
    articulation_points,

    Chain,
    run_chain,

    # swap up/swap down
    pop_balance, 
    all_dist_art_pts, 
    project_up_pctpair, 
    project_down_pctpair, 
    count_split_nodes, 
    get_valid_splits_pair_tempered, 
    get_valid_merges_pair, 
    iso_if_you_flipped, 
    swap_pair_pct_tempered, 
    copy_partition, 
    swap_up, 
    swap_down, 

    # tempering hierarchy 
    TemperingHierarchy, 
    add_level!, 
    write_TH_partition_to_file,
    highest_ancestor,
    ancestor_at_level,
    delete_level_at_top!,
    PT_swap!,

    generate_level, 

    # random utils 
    reverse_vmap,
    uniformly_choose,

    log_nspanning,

    populate_extras!,
    remove_extras!,
    check_onbr_numbers_match,
    write_partition_to_file_flat,
    remove_artpt_forming_edges

include("./graph.jl")
include("./SimpleWeightedGraphs_BugFixes.jl")
include("./constraint_types.jl")

include("./utils.jl")
include("./multi_level_graph.jl")
include("./expand_hierarchy.jl")
include("./multi_level_partition.jl")

include("./energy_types.jl")
include("./measure.jl")

include("./simple_writer.jl")


include("./cuttable_tree.jl")
include("./forest_recom2.jl")
include("./multistep_forest_recom2.jl")
include("./single_node_flip.jl")
include("./mcmc.jl")
include("./constraints.jl")

include("./multistep_SNF.jl")
include("./chain.jl")
include("./parallel_tempering.jl")

include("./polsby_popper.jl")
include("./vap_frac.jl")

include("quad_edge.jl")
include("hierarchical_partition.jl")

include("./tempering_hierarchy.jl")
include("./hierarchy_swapping.jl")
include("./hierarchy_gen.jl")

# include("./parallel_tempering_multiprocessing.jl")
include("./tree.jl")

end # module
