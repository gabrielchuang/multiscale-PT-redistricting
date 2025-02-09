using StructTypes,JSON3,Dates 

struct SimpleWriter
    io::IO
end 

function SimpleWriter(
    measure::Measure,
    constraints::Dict,
    partition::MultiLevelPartition,
    output_file_path::String;
    description::String="",
    time_stamp=string(Dates.now()),
    rng_label::String=""
)
    SimpleWriter(measure, constraints, partition.num_dists, output_file_path; 
        description=description, time_stamp=time_stamp, rng_label=rng_label)
end

function SimpleWriter(
    measure::Measure,
    constraints::Dict,
    num_dists::Int,
    output_file_path::String;
    description::String="",
    time_stamp=string(Dates.now()),
    rng_label::String=""
)
    parameters = Dict(
        "gamma"=>measure.gamma, 
        "energies"=>measure.descriptions,
        "energy weights"=>measure.weights,
        "districts"=>num_dists
    )
    if haskey(constraints, PopulationConstraint)
        min_pop = constraints[PopulationConstraint].min_pop
        max_pop = constraints[PopulationConstraint].max_pop
        parameters["population bounds"] = [min_pop, max_pop]
    end

    header = Dict("date" => time_stamp, description => description)

    dir = dirname(output_file_path)
    if !isdir(dir)
        mkpath(dir)
    end

    io = smartOpen(output_file_path, "w")

    JSON3.write(io, "redistricting maps: "*rng_label)
    write(io,"\n")
    JSON3.write(io,header)
    write(io,"\n")
    JSON3.write(io,parameters)
    write(io,"\n")

    return SimpleWriter(io)
end

function write_districts(
    writer::SimpleWriter, 
    partition::MultiLevelPartition, 
    d_initial::Dict = Dict(); 
    use_legacy_output=true
)
    use_legacy_output = false
    if use_legacy_output
        d = [] #Dict{Vector{String}, Int}()
        for (k,v) in partition.node_to_district 
            push!(d, 
            Dict([node_name_exclude_quotient_nodes(partition.graph, k) => v])
            )
        end
        d_initial["districting"] = d 
        buff=JSON3.write(d_initial)
        write(writer.io, buff)
        write(writer.io, "\n") 

    else 
        d = Dict{GlobalID, Vector{GlobalID}}() 
        for i = 1:partition.num_dists
            d[i] = [k for (k,v) in partition.district_to_nodes[i] if v === nothing]
        end 
        d_initial["districts"] = d 
        buff=JSON3.write(d_initial)
        write(writer.io, buff)
        write(writer.io, "\n")
    end 
end 

function close_writer(writer::SimpleWriter)
    close(writer.io)
end


function smartOpen(fileName::String, io_mode::String)::Union{IO,Nothing}
    ext,base=getFileExtension(fileName)
 
    if ((io_mode=="w") |(io_mode=="a"))
        
        oo= try open(fileName,io_mode) #try to open filename given
        catch err
            @info( string("Error opening file. Trying alternative extensions for ",fileName));
            if ((io_mode=="a") & ((ext==".gz") | (ext==".bz2")))
                fileName=base;
                ext,base=getFileExtension(base);
                open(base,io_mode) #if first open fails and was compressed, try not compressed
            else
                base=string(base,ext);
                ext=".gz";
                fileName=string(base,ext);
                open(fileName,io_mode) #if first open fails and was not compressed, try  compressed
            end
        end
        if ext==".bz2"
            # print("w-bZ")
            return Bzip2CompressorStream(oo)
        end
        if ext==".gz"
            # print("w-gZ")
            return GzipCompressorStream(oo)
        end
        return oo  
    end
    
    if io_mode=="r" 
        
            oo= try open(fileName,io_mode) #try to open filename given
             catch err
                @info( string("Error opening file. Trying alternative extensions for ",fileName));
                if ((ext==".gz") | (ext==".bz2"))
                    fileName=base;
                    ext,base=getFileExtension(base);
                    oo=open(fileName,io_mode) #if first open fails and was compressed, try not compressed
                else
                    base=string(base,ext);
                    ext=".gz";
                    fileName=string(base,ext);
                    open(fileName,io_mode) #if first open fails and was not compressed, try  compressed
            end
        end
        if ext==".bz2"
            # print("r-bZ")
            return Bzip2DecompressorStream(oo)
        end
        if ext==".gz"
            # print("r-gZ")
            return GzipDecompressorStream(oo)
        end
        return oo
    end
    
    if ext==".bz2"
            return nothing
    end
    if ext==".gz"
            return nothing
    end
    oo=open(fileName,io_mode)
    return oo
end
function getFileExtension(filename::String)
    i=findlast(isequal('.'),filename)
    return filename[i:end],filename[1:i-1]
end

