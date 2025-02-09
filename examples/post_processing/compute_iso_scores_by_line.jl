import Pkg
Pkg.activate("mergeSplit"; shared=true)

import Base.Threads.@spawn
using Dates, JSON
using MultiScaleMapSampler
using CSV,Tables,LazyArrays,Base#,Revise

print("Starting processing (",Dates.now(),")...\n");flush(stdout) 

# expected directory structure: same as compute_marginals_by_line

function get_garea(g, nodes)
    return sum([get_node_area(g, x) for x in nodes])
end

function get_gperimeter(g, nodes)
    total_perim = 0 
    nodes = Set(nodes)
	for n in nodes 
        total_perim += node_attributes(g,n)["border_length"]
        for (n2, attrs) in get_finest_level_neighbors(g, n)
            if n2 in nodes
                continue 
            else 
                total_perim += attrs.length 
            end
        end
    end
    return total_perim 
end 

identifier = ARGS[1]

for identifier in ARGS

    input_dir = joinpath("..", "output", identifier*"_districts")
    output_dir=joinpath("..", "output", identifier*"_marginals") 
    #pct_graph_path = joinpath("..", "test_graphs", "NC_pct21.json")
    #pct_graph_path = joinpath("..", "test_graphs", "IA_pct20.json")
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
    nodeData = Set(["COUNTY", "NAME", "area", "border_length", "POP20", "G20PREDEM", "G20PREREP"]);
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


        io = open(joinpath(input_dir, filename), "r")

        ct = 0
        buff=readline(io) 
        iso_scores = [] 
        while buff != "" 
            districting = JSON.parse(buff)
            if typeof(districting) == String || !haskey(districting, "name")
                # invalid row (parameters, etc)
            else 
                this_plan_scores = [] 
                for (_, leaf_nodes) in districting["districts"]
                    dist_area = get_garea(graph, leaf_nodes) 
                    dist_perimeter = get_gperimeter(graph, leaf_nodes) 
                    push!(this_plan_scores, dist_perimeter^2 / dist_area)
                end 
                sort!(this_plan_scores)
                push!(iso_scores, this_plan_scores)
            end 

            buff=readline(io) 
            ct += 1 

            if mod(ct, 10000) == 0
                @show ct 
            end 
        end

        #tc_file = smartOpen(joinpath(output_dir_prefix, file_desc*"_forest_counts.csv"), "w")
        #write(tc_file, join(forest_counts, ", "))

        iso_file = open(joinpath(output_dir_prefix, file_desc*"_iso_scores.csv"), "w")
        write(iso_file, "iso_scores\n"*join([join(x, ",") for x in iso_scores], "\n"))

        @show ct 
    end
end 
