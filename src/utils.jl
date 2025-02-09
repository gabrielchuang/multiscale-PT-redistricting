toggle_ENSURE() = true # when set to true, this will enable all `@ensure`s

macro ensure(test)
  esc(:(if $(@__MODULE__).toggle_ENSURE()
    @assert($test)
   end))
end

function reverse_vmap(vmap) 
    Dict([vmap[x] => x for x in keys(vmap)])
end

function select(p, collection) 
    @ensure length(collection) > 0 
    return first(x for x in collection if p(x)) 
end 

function uniformly_choose(
    set_to_choose_from,
    rng::AbstractRNG
)
    if length(set_to_choose_from)==0
        return nothing, nothing
    end
    rnd_ind = Int(ceil(rand(rng)*length(set_to_choose_from)))
    return collect(set_to_choose_from)[rnd_ind], 1.0/length(set_to_choose_from)
end

function max_index(coll, sz) 
    max_sz = 0
    max_idx = 0
    for i in eachindex(coll)
        if sz(coll[i]) > max_sz 
            max_sz = sz(coll[i])
            max_idx = i 
        end
    end
    return max_idx 
end

function set_matrix!(matrix, contents) 
    @assert size(matrix) == size(contents) 
    n,m = size(matrix)
    for i=1:n 
        for j=1:m 
            matrix[i,j] = contents[i,j]
        end
    end
end

# to keep structs immutable, we sometimes have to do this to their fields 
# (rather than struct.coll = new_contents, which isn't allowed 
# for immutable structs). 
function replace!(coll, new_contents) 
    empty!(coll)
    append!(coll, new_contents) 
end 

function print_dict(dict::Dict, line_prefix="")
    for (k,v) in dict 
        if typeof(v) <: Dict 
            println(line_prefix, k)
            print_dict(v, line_prefix*" ")
        elseif typeof(k) <: Dict 
            
        elseif length(v) > 10 
            println(line_prefix, k, ": ", "...")
        else
            println(line_prefix, k, ": ", v)
        end
    end 
end 
