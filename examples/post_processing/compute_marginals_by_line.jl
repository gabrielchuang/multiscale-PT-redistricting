import Pkg
Pkg.activate("mergeSplit"; shared=true)

import Base.Threads.@spawn
using Dates, JSON
using MultiScaleMapSampler
using CSV,Tables,LazyArrays,Base#,Revise

print("Starting processing (",Dates.now(),")...\n");flush(stdout) 

#= expected directory structure: 
test/ 
    output/ 
        identifier_districts/ 
            identifier_hierarchy.json # unclear 
            each_run1.json 
            each_run2.json 
            ...
        identifier_marginals/ (creates) 
            each_run1_marginals/ 
                each_run1_marginals_election1.csv
                each_run2_marginals_election1.csv
                ...
            ... 
    test_graphs/NC_pct21.json 
    post_processing/this_file.jl 
=#

# tail -q -n +4 $(ls -t NC_cluster_500k_A*.jsonl) > 500k_A.jsonl
# tail -q -n +4 $(ls -t CT4_10_E*_2.jsonl) > CT4_10_E.jsonl

identifier = ARGS[1]

for identifier in ARGS

    input_dir = joinpath("..", "output", identifier*"_districts")
    output_dir=joinpath("..", "output", identifier*"_marginals") 
    #pct_graph_path = joinpath("..", "test_graphs", "NC_pct21.json")
    pct_graph_path = joinpath("..", "test_graphs", "CT_pct20.json")
    #pct_graph_path = joinpath("..", "test_graphs", "coarse_graphs_comp", "NC_16d_blocks.json")
    #hierarchy_path = joinpath(input_dir, identifier*"_hierarchy.json")

    print("Reading from:\t", input_dir,"\n")
    print("Writing to:\t", output_dir,"\n"); flush(stdout)

    files = [f for f in readdir(input_dir) if occursin("jsonl", f)]

    #elections=["G16_PR"]#, "G16_USS", "G16_GV", "G16_LG", "G16_AG"]; 
    #variables=["pop2020cen"];

    #pop_field = "pop2020cen"
    #nodeData = Set(["county", "prec_id", "area", pop_field]);
    nodeData = Set(["COUNTY", "NAME", "area", "POP20", "G20PREDEM", "G20PREREP"]);
    #nodeData = Set(["block_name", "area", "id", pop_field]);
    #nodeData = union(nodeData, Set([string(e,"_R") for e in elections ]))  #  load election republican votes
    #nodeData = union(nodeData, Set([string(e,"_D") for e in elections ]))  #  load election democratic votes
    #nodeData = union(nodeData, Set(variables))                             #  load additional variables

    base_graph = BaseGraph(pct_graph_path, "POP20", inc_node_data=nodeData, area_col="area")
    for k in base_graph.node_attributes 
        k["county_and_pct"] = string(k["COUNTY"], ",",k["NAME"])
    end
    

    #graph = MultiLevelGraph(base_graph, ["block_name"]); 
    #graph = MultiLevelGraph(base_graph, ["prec_id", "county"], hierarchy_path); 
    #graph = MultiLevelGraph(base_graph, ["county_and_pct"], hierarchy_path); 
    graph = MultiLevelGraph(base_graph, ["county_and_pct"]); 

    print("Will process: \n")
    for a in files
        println("\t",a)
    end

    for filename in files
        println("processing: ", filename)

        file_desc = split(filename, ".")[1]
        output_dir_prefix=joinpath(output_dir,file_desc*"_marginals")
        mkpath(output_dir_prefix)

        output_files_by_election = Dict()
        for election in elections
            out_io = open(joinpath(output_dir_prefix, file_desc*"_marginals_"*election*".csv"), "w")
            output_files_by_election[election] = out_io
        end 

        io = open(joinpath(input_dir, filename), "r")

        ct = 0
        buff=readline(io) 
        forest_counts = [] 
        iso_scores = [] 
        while buff != "" 
            districting = JSON.parse(buff)
            if !haskey(districting, "name")
                # invalid row (parameters, etc)
            else 
                for election in elections 
                    district_vote_shares = []
                    for (di, leaf_nodes) in districting["districts"] 
                        votes_D = sum([node_attributes(graph, node_ID)[election*"_D"] 
                            for node_ID in leaf_nodes])
                        votes_R = sum([node_attributes(graph, node_ID)[election*"_R"] 
                            for node_ID in leaf_nodes])
                        
                        district_vote_share = votes_D / (votes_D + votes_R)
                        push!(district_vote_shares, district_vote_share)
                    end
                    dvs_string = join(district_vote_shares, ", ")
                    write(output_files_by_election[election], string(districting["step"])*",")
                    write(output_files_by_election[election], dvs_string)
                    write(output_files_by_election[election], "\n")
                end

                # the below only works if it is a pct-only graph 
                #=log_forest_count = 0 
                iso_score = 0 
                for (_, leaf_nodes) in districting["districts"]
                    mask = Dict{Int64, Union{Nothing, Set{Int64}}}([x => nothing for x in leaf_nodes])
                    mask[graph.root.global_ID] = Set{Int64}(leaf_nodes) 

                    dist_subgraph = take_subgraph(graph, mask)
                    #precompute_log_tree_counts!(dist_subgraph) 
                    #log_forest_count += log_hierarchical_tree_count(dist_subgraph)

                    dist_area = get_area(dist_subgraph) 
                    dist_perimeter = get_perimeter(graph, leaf_nodes) 
                    iso_score += dist_perimeter^2 / dist_area 
                end 
                #push!(forest_counts, log_forest_count)
                push!(iso_scores, iso_score)=#
            end 

            buff=readline(io) 
            ct += 1 

            if mod(ct, 10000) == 0
                @show ct 
            end 
        end

        #tc_file = smartOpen(joinpath(output_dir_prefix, file_desc*"_forest_counts.csv"), "w")
        #write(tc_file, join(forest_counts, ", "))

        #iso_file = smartOpen(joinpath(output_dir_prefix, file_desc*"_iso_scores.csv"), "w")
        #write(iso_file, "iso_scores\n"*join(iso_scores, "\n"))

        @show ct 
    end
end 
