using MultiScaleMapSampler
using Test
using RandomNumbers

const testdir = dirname(@__FILE__)

small_square_json = joinpath("test_graphs", "4x4pct_2x2cnty.json")
small_square_node_data = Set(["county", "pct", "pop", "area", "border_length"])
6small_square_base_graph = BaseGraph(small_square_json, "pop", inc_node_data=small_square_node_data,
                                    area_col="area",node_border_col="border_length", 
                                    edge_perimeter_col="length")
small_square_graph = MultiLevelGraph(small_square_base_graph, ["pct", "county"])
small_square_dists = 2

three_level_json = joinpath("test_graphs", "2x6pct_3level.json")
three_level_node_data = Set(["county", "pct", "cblock", "pop", "area", "border_length"])
three_level_base_graph = BaseGraph(three_level_json, "pop", inc_node_data=three_level_node_data,
                                    area_col="area",node_border_col="border_length", 
                                    edge_perimeter_col="length")
three_level_graph = MultiLevelGraph(three_level_base_graph, ["cblock", "pct", "county"])
three_level_dists = 2

tests = [
    "small_square_p88",
    "small_square_p79", 
    "small_square_p79_comp"
    ]

for t in tests
    tp = joinpath(testdir, "test_cases/$(t).jl")
    include(tp)
end