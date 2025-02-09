mutable struct Chain{T <: Real}
    proposal::Union{Function,Vector{Tuple{T, Function}}}
    measure::Measure
    writer::Union{SimpleWriter, Nothing}
    rng::AbstractRNG
end

function Chain(
	proposal::Function,
	measure::Measure, 
	writer::Union{SimpleWriter, Nothing},
    rng::AbstractRNG
)
	return Chain{Real}(proposal, measure, writer, rng)
end

function run_chain!(
	partition::MultiLevelPartition, 
	chain::Chain,
	steps::Union{Int,Tuple{Int,Int}}
)
	try
		run_metropolis_hastings!(partition, chain.proposal, chain.measure, steps, 
			                     chain.rng, writer=chain.writer)
		return 0, chain.measure.gamma
	catch e
		return 1, chain.measure.gamma
	end
end