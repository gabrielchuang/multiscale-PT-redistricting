import Base.Threads.@spawn
# note: this file has not been tested or used yet, merely 
# refactored to use the new multi level graph/partition structure. 

function parallel_tempering_multithread!(
	replicas::Vector{MultiLevelPartition},
	chains::Vector{Chain},
	steps::Int,
	swap_interval::Int
)
	threads = Threads.nthreads()
	@assert length(replicas) == length(chains)
	@assert length(replicas) == threads
	initial_step, final_step = 1, steps

	swaps = floor((final_step-initial_step)/swap_interval)

	# extend partitions
	for ii = 1:length(replicas)
		replicas[ii].extensions[replica_id::EXTENSIONS] = ii
		replicas[ii].extensions[bath_swaps::EXTENSIONS] = 0
		num_dists = replicas[ii].num_dists
		replicas[ii].extensions[del_dists::EXTENSIONS]=zeros(num_dists)
	end
	error_msg = zeros(Threads.nthreads(), 2)

	pair_start, no_tasks = 0, 0
	tasks = Vector(undef, length(replicas))
	# @show length(tasks)

	for swap = 1:swaps
		tasks = []
		for ii = 1:length(replicas)
			inc = (Int((swap-1)*swap_interval), Int(swap*swap_interval))
			push!(tasks, @spawn run_chain!(replicas[ii], chains[ii], inc))
		end
		for ii = 1:length(replicas)
		 	msg = fetch(tasks[ii])
		 	error_msg[ii, 1] = msg[1]
		 	error_msg[ii, 2] = msg[2]
		end
		if maximum(error_msg[:, 1]) > 0
			println("erroring out")
			@show error_msg
			break
		end

		println("swap ", swap, " out of ", swaps); flush(stdout)
		pair_start += 1
		task_ind = 1
		tasks = []
		
		for ii = pair_start:2:length(replicas)-1
			# tasks[task_ind] = Base.Threads.@spawn try_swap_replicas!(replicas, chains, ii)
			push!(tasks, @spawn try_swap_replicas!(replicas, chains, ii))
			task_ind += 1
		end
		for t in tasks
		 	wait(t)
		end

		pair_start = mod(pair_start, 2)
	end
end

""""""
function try_swap_replicas!(
	replicas::Vector{MultiLevelPartition},
	chains::Vector{Chain},
	index::Int
)::Nothing 
	@assert index < length(chains)
	@assert length(replicas) == length(chains)

	log_p_ii = log_measure(replicas[index], chains[index].measure)
	log_p_jj = log_measure(replicas[index+1], chains[index+1].measure)
	log_p_ij = log_measure(replicas[index], chains[index+1].measure)
	log_p_ji = log_measure(replicas[index+1], chains[index].measure)
	accept_prob = exp(log_p_ji + log_p_ij - log_p_ii - log_p_jj)

	println("try_swap_replicas: ", index, " ", accept_prob, " ", [log_p_ii, log_p_jj, log_p_ij, log_p_ji])
	if rand(chains[index].rng) < accept_prob
		tmp = replicas[index]
		replicas[index] = replicas[index+1]
		replicas[index+1] = tmp
	end

	# ??????? these are non-state-affecting functions. what is the purpose of calling them here? 
	get_log_energy(replicas[index], chains[index].measure) 
	get_log_energy(replicas[index+1], chains[index+1].measure)
end

function parallel_tempering_multiprocess!(
	replicas::Vector{MultiLevelPartition},
	chains::Vector{Chain},
	steps::Int,
	swap_interval::Int
)
end 

#using MPI

function parallel_tempering!(
    replica::MultiLevelPartition,
    chain::Chain,
    steps::Int,
    swap_interval::Int,
    measure_parameters::Array{Float64},
    measure_to_rank::Vector{Int},
    proposal_weights::Array{Float64},
    base_sampler::Union{Nothing, Tuple{Measure, Vector}}=nothing
)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    tasks = MPI.Comm_size(comm)

    initial_step, final_step = set_step_bounds(steps)
    swaps = Int(floor((final_step-initial_step)/swap_interval))

#######
    # extend partitions
    num_dists = replica.num_dists
    replica.extensions[del_dists::EXTENSIONS] = zeros(num_dists)
    replica.extensions[bath_swaps::EXTENSIONS] = 0
    replica.extensions[measure_id::EXTENSIONS] = rank+1
#######
    pair_start = 0

    for swap = 1:swaps
        inc = (Int((swap-1)*swap_interval), Int(swap*swap_interval))
        t1 = time_ns()
        error_msg = run_chain!(replica, chain, inc)
        t2 = time_ns()
        if error_msg[1] > 0
            println("erroring out")
            @show error_msg
            @assert false
        end

        for jj = 0:tasks-1
            if jj == rank
                if rank == 0
                    println("swap ", swap, " out of ", swaps); flush(stdout)
                end
                m_id = replica.extensions[measure_id::EXTENSIONS]
                time_run = (t2-t1)/10^9
                @show rank, m_id, time_run
                flush(stdout)
            end
            MPI.Barrier(comm)
        end
        MPI.Barrier(comm)
        swap_energies!(measure_to_rank, replica, chain, pair_start, swap, 
                       measure_parameters, proposal_weights, base_sampler)
        flush(stdout); MPI.Barrier(comm)

        pair_start = Int(!Bool(pair_start)) # 0 -> 1 -> 0
    end
end

