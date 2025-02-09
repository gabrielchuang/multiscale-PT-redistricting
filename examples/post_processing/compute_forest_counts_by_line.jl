import Pkg
Pkg.activate("mergeSplit"; shared=true)

#import Base.Threads.@spawn
using Dates, JSON
using MultiScaleMapSampler
using CSV,Tables,Base#,Revise

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
    output_dir=joinpath("..", "output", identifier*"_observables") 
    pct_graph_path = joinpath("..", "test_graphs", "16x16pct.json")

    print("Reading from:\t", input_dir,"\n")
    print("Writing to:\t", output_dir,"\n"); flush(stdout)

    files = [f for f in readdir(input_dir) if occursin("jsonl", f)]

    nodeData = Set(["county", "pct", "pop", "area", "border_length", "id"])

    base_graph = BaseGraph(pct_graph_path, "pop", inc_node_data=nodeData, area_col="area")
    for k in base_graph.node_attributes 
        k["county_and_pct"] = string(k["county"], ",",k["pct"])
    end
    
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

        io = open(joinpath(input_dir, filename), "r")

        ct = 0
        buff=readline(io) 
        forest_counts = [] 
        district_tree_counts = []
        iso_scores = [] 
        while buff != "" 
            districting = JSON.parse(buff)
            if typeof(districting) >: String || !haskey(districting, "name")
                # invalid row (parameters, etc)
            else 
                # the below only works if it is a pct-only graph 
                tree_counts = [] 
                log_forest_count = 0 
                for (_, leaf_nodes) in districting["districts"]
                    #@show leaf_nodes
                    mask = Dict{Int64, Union{Nothing, Set{Int64}}}([x => nothing for x in leaf_nodes])
                    mask[graph.root.global_ID] = Set{Int64}(leaf_nodes) 
                    dist_subgraph = take_subgraph(graph, mask)
                    
                    s_graph = child_graph(dist_subgraph, dist_subgraph.root.global_ID)[1]
                    lnsp = log_nspanning(s_graph)
                    #@show lnsp 
                    
                    precompute_log_tree_counts!(dist_subgraph) 


                    #@show dist_subgraph.log_tree_counts
                    #@show dist_subgraph.root.global_ID 
                    lhtc = log_hierarchical_tree_count(dist_subgraph)
                    log_forest_count += lhtc
                    push!(tree_counts, lhtc)
                    #@show lhtc
                end 
                push!(forest_counts, log_forest_count)
                push!(district_tree_counts, join(tree_counts, ","))
            end 

            buff=readline(io) 
            ct += 1 

            if mod(ct, 10000) == 0
                @show ct 
            end 
        end

        fc_file = open(joinpath(output_dir_prefix, file_desc*"_forest_counts.csv"), "w")
        write(fc_file, join(forest_counts, ", "))

        tc_file = open(joinpath(output_dir_prefix, file_desc*"_tree_counts.csv"), "w")
        write(tc_file, join(district_tree_counts, "\n"))

		@show ct 
    end
end 
