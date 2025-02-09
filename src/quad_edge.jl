const QuadEdgeID = Int 
const QuadFaceID = Int 


# the dual of the mesh is obtained by swapping vertices with faces, 
# src with left_face, and left_face with ccw_edge.
const MultiQuadEdge = Union{QuadEdgeID, Set{QuadEdgeID}}

struct QuadEdgeMesh 
    # for each vertex v, we hold one of the edges e with e.src = v 
    vertices::Dict{GlobalID, QuadEdgeID} 
    # for each face f, we hold one of the edges with e.left_face = f 
    faces::Dict{QuadFaceID, QuadEdgeID}
    # given src, dst, find the quadedge 
    edge_lookup::Dict{Tuple{GlobalID, GlobalID}, MultiQuadEdge}

    # finally, each quadedge holds a few things: 
    srcs::Dict{QuadEdgeID, GlobalID} # the src of the edge 
    syms::Dict{QuadEdgeID, QuadEdgeID} # the reverse edge 
    left_face::Dict{QuadEdgeID, QuadFaceID} # left face
    onext::Dict{QuadEdgeID, QuadEdgeID} # the next edge, clockwise, around the same src    
end

function QuadEdgeMesh(
    G::MultiLevelGraph, 
    oriented_neighbor_list::Dict{GlobalID, Vector{GlobalID}} # clockwise list 
)
    return QuadEdgeMesh(finest_level_nodes(G), oriented_neighbor_list)
end

function next_element_in_nodes(onbrs, nodes, i) 
    while true #c < 100 
        i += 1 
        if i > length(onbrs) 
            i = 1 
        end 
        if onbrs[i] in nodes
            return onbrs[i] 
        end
    end
end 

function QuadEdgeMesh(
    nodes::Set{GlobalID}, 
    oriented_neighbor_list::Dict{GlobalID, Vector{GlobalID}} # clockwise list 
)
    # initialize the edge IDs and onext 
    edge_lookup = Dict{Tuple{GlobalID, GlobalID}, Union{QuadEdgeID, Set{QuadEdgeID}}}()
    vertices = Dict{GlobalID, QuadEdgeID}()
    srcs = Dict{QuadEdgeID, GlobalID}()
    dsts = Dict{QuadEdgeID, GlobalID}()
    i = 0 
    onext = Dict{QuadEdgeID, QuadEdgeID}()
    oprev = Dict{QuadEdgeID, QuadEdgeID}()
    # invariant: onext[oprev[x]] == x and oprev[onext[x]] == x. 
    for src in nodes
        cur_i = i

        if length(oriented_neighbor_list[src]) == 0
            println("node has no oriented nbrs: ", src) 
        end

        zero_onbrs = true 
        for dst in oriented_neighbor_list[src]
            if !(dst in nodes) continue end 
            if haskey(edge_lookup, (src, dst))
                #println("two pcts border in more than one place") 
                if typeof(edge_lookup[(src,dst)]) >: QuadEdgeID
                    edge_lookup[(src,dst)] = Set([edge_lookup[(src,dst)], i])
                else 
                    push!(edge_lookup[(src,dst)], i)
                end
            else
                edge_lookup[(src, dst)] = i
            end
            onext[i] = i+1
            oprev[i] = i-1
            srcs[i] = src 
            dsts[i] = dst
            vertices[src] = i 
            i += 1
            zero_onbrs = false 
        end
        @assert !zero_onbrs
        onext[i-1] = cur_i
        oprev[cur_i] = i-1
    end

    for (k,v) in onext 
        if oprev[v] != k || onext[k] != v
            @show oprev[v], k, onext[k], v
            @show onext[v], oprev[k]
            @show onext[onext[v]], oprev[onext[v]]
            @show oprev[oprev[k]], onext[oprev[k]]
        end
        @assert oprev[v] == k 
        @assert onext[k] == v
    end 

    sym = Dict{QuadEdgeID, QuadEdgeID}() 
    # set sym 
    for src in nodes
        for dst in oriented_neighbor_list[src]
            if !(dst in nodes) continue end 
            edge = edge_lookup[(src, dst)] 
            if haskey(sym, edge) continue end 

            if typeof(edge) >: QuadEdgeID
                if !(typeof(edge_lookup[(dst, src)]) >: QuadEdgeID)
                    @show dst, src 
                    @assert false 
                end
                sym[edge] = edge_lookup[(dst, src)]
                sym[edge_lookup[(dst, src)]] = edge 
            else 
                # this is a multiedge. how to tell which is its sym?? :( 
                @assert typeof(edge) >: Set{QuadEdgeID}
                if !(typeof(edge_lookup[(dst, src)]) >: Set{QuadEdgeID})
                    @show typeof(edge_lookup[(dst, src)])
                    @show dst, src 
                    @assert false "sym asymmetric" 
                end
                multiedge = edge
                
                # find a common neighbor that both src and dst touch only once. 
                # then, going clockwise for one and ccw for the other, enumerate in order 
                # and pair up the ones that come in sequence. 
                single_touch_srcs = filter(x -> count(y -> y == x, oriented_neighbor_list[src]) == 1 && x in nodes, oriented_neighbor_list[src])
                single_touch_dsts = filter(x -> count(y -> y == x, oriented_neighbor_list[dst]) == 1 && x in nodes, oriented_neighbor_list[dst])

                src_to_common = nothing; dst_to_common = nothing 

                if length(intersect(Set(single_touch_srcs), Set(single_touch_dsts))) == 0 
                    @show dst, src, single_touch_srcs, single_touch_dsts 
                    #bleh. we'll use multiple touch points and hope that planarity keeps things sane 
                    single_touch_srcs = filter(x -> x in nodes, oriented_neighbor_list[src])
                    single_touch_dsts = filter(x -> x in nodes, oriented_neighbor_list[dst])    

                    common_single_touch = first(intersect(Set(single_touch_srcs), Set(single_touch_dsts)))
                    dst_to_common = edge_lookup[(dst, common_single_touch)] 
                    src_to_common = edge_lookup[(src, common_single_touch)]
                    
                    if typeof(src_to_common) >: Set{QuadEdgeID} 
                        src_to_common = first(src_to_common) 
                    end 
                    if typeof(dst_to_common) >: Set{QuadEdgeID} 
                        dst_to_common = first(dst_to_common)
                    end
                else 
                    common_single_touch = first(intersect(Set(single_touch_srcs), Set(single_touch_dsts)))
                    dst_to_common = edge_lookup[(dst, common_single_touch)] 
                    src_to_common = edge_lookup[(src, common_single_touch)]
                end
                @assert typeof(dst_to_common) >: QuadEdgeID 
                @assert typeof(src_to_common) >: QuadEdgeID 

                for _ = 1:length(multiedge)
                    # go clockwise from dst_to_common and ccw from src_to_common

                    while dsts[dst_to_common] != src 
                        dst_to_common = oprev[dst_to_common]
                    end 
                    while dsts[src_to_common] != dst 
                        src_to_common = onext[src_to_common]
                    end 
                    sym[dst_to_common] = src_to_common
                    sym[src_to_common] = dst_to_common
                    dst_to_common = oprev[dst_to_common]
                    src_to_common = onext[src_to_common]
                end
            end
        end
    end

    #println("syms have been set")
    #@assert false "spot 1"

    for (_, edge) in edge_lookup 
        if typeof(edge) >: Set{QuadEdgeID}
            for e in edge 
                if sym[sym[e]] != e 
                    @show e 
                    @assert false 
                end
                @assert oprev[onext[e]] == e 
                @assert onext[oprev[e]] == e
            end 
        else 
            @assert sym[sym[edge]] == edge 
            @assert oprev[onext[edge]] == edge 
            @assert onext[oprev[edge]] == edge
        end
    end

    #println("so far, all checks ok")
    #@assert false "spot 2"
   
    # set lface and faces 
    lface = Dict{QuadEdgeID, QuadFaceID}()
    faces = Dict{QuadFaceID, QuadEdgeID}()
    i = 0
    for src in nodes
        for dst in oriented_neighbor_list[src]
            if !(dst in nodes) continue end 
            cur_edge = edge_lookup[(src, dst)]
            if typeof(cur_edge) >: Set{QuadEdgeID}
                #cur_edge = first(cur_edge) 

                for ce2 in cur_edge
                    if haskey(lface, ce2) 
                        continue 
                    end 
        
                    left_trav = traverse_face(onext, sym, ce2)
                    for lt in left_trav 
                        lface[lt] = i 
                    end 
                    faces[i] = ce2 
                    i += 1 
                end
            else
                if haskey(lface, cur_edge) 
                    continue 
                end 

                left_trav = traverse_face(onext, sym, cur_edge)
                for lt in left_trav 
                    lface[lt] = i 
                end 
                faces[i] = cur_edge 
                i += 1 
            end
        end
    end

    n_edges = sum([typeof(x) >: Set{QuadEdgeID} ? length(x) : 1 for (_, x) in edge_lookup])
    if length(faces) + length(vertices) - n_edges/2 != 2.0
        println("multiple connected components in this QEM: ", (length(faces) + length(vertices) - n_edges/2)/2)
    end
    return QuadEdgeMesh(vertices, faces, edge_lookup, srcs, sym, lface, onext)
end


# return all the edges, in order, that have the same left face as cur_edge
function traverse_face(onext, sym, cur_edge) 
    left_trav = [] 
    og_cur_edge = cur_edge 
    for emergency_brake = 1:1000
        push!(left_trav, cur_edge)
        cur_edge = onext[sym[cur_edge]]
        if (cur_edge == og_cur_edge) 
            return left_trav 
        end
    end
    @assert false "somehow looping forever; sym or onext set incorrectly"
end

function traverse_vertex(onext, cur_edge)
    v_trav = [] 
    og_cur_edge = cur_edge 
    for emergency_brake = 1:1000
        push!(v_trav, cur_edge)
        cur_edge = onext[cur_edge]
        if (cur_edge == og_cur_edge) 
            return v_trav 
        end
    end
    @assert false "somehow looping forever; sym or onext set incorrectly"
end

function delete_edge!(QEM, e1)
    e2 = QEM.syms[e1]

    f1 = QEM.left_face[e1]
    f2 = QEM.left_face[e2]
    left_trav_f1 = traverse_face(QEM.onext, QEM.syms, e1)
    left_trav_f2 = traverse_face(QEM.onext, QEM.syms, e2)
    v_trav_e1 = traverse_vertex(QEM.onext, e1)
    v_trav_e2 = traverse_vertex(QEM.onext, e2)

    for left_face_edge in left_trav_f1 
        QEM.left_face[left_face_edge] = f2 
    end

    e1_onext = QEM.onext[e1] 
    e1_oprev = v_trav_e1[end]
    e2_onext = QEM.onext[e2]
    e2_oprev = v_trav_e2[end]
    QEM.onext[e1_oprev] = e1_onext
    QEM.onext[e2_oprev] = e2_onext


    if QEM.faces[QEM.left_face[e1]] == e1 
        QEM.faces[QEM.left_face[e1]] = left_trav_f1[end]
    end 
    if QEM.faces[QEM.left_face[e2]] == e2 
        QEM.faces[QEM.left_face[e2]] = left_trav_f2[end] 
    end 

    if QEM.vertices[QEM.srcs[e1]] == e1 
        QEM.vertices[QEM.srcs[e1]] = e1_onext
    end 
    if QEM.vertices[QEM.srcs[e2]] == e2
        QEM.vertices[QEM.srcs[e2]] = e2_onext
    end 

    delete!(QEM.edge_lookup[(QEM.srcs[e1], QEM.srcs[e2])], e1)
    delete!(QEM.edge_lookup[(QEM.srcs[e2], QEM.srcs[e1])], e2)

    delete!(QEM.onext, e1)
    delete!(QEM.onext, e2)
    delete!(QEM.left_face, e1)
    delete!(QEM.left_face, e2)
    delete!(QEM.srcs, e1)
    delete!(QEM.srcs, e2)
    delete!(QEM.syms, e1)
    delete!(QEM.syms, e2)

    delete!(QEM.faces, f1)
end

function remove_multiedges!(QEM::QuadEdgeMesh) 
    for (k,v) in QEM.edge_lookup
        if typeof(v) >: Set{QuadEdgeID} 
            # this is a multi-edge. remove all but one of them and relabel the faces 
            # (as well as onext/oprev)
            v2 = collect(v)
            for edge in v2[2:end]
                # replace all instances of f1 with f2 
                #println("removing edge ", edge, " = ", QEM.srcs[edge], " - ", QEM.srcs[QEM.syms[edge]])
                delete_edge!(QEM, edge)
            end
            QEM.edge_lookup[k] = v2[1]
        end
    end
end

function edge_lookup_by_face(QEM::QuadEdgeMesh, f1::QuadFaceID, f2::QuadFaceID)::QuadEdgeID
    # returns the edge E such that E.left_face = f1 and E.sym.left_face = f2. 
    cur_edge = QEM.faces[f1]
    while QEM.left_face[QEM.syms[cur_edge]] != f2 
        cur_edge = QEM.onext[QEM.syms[cur_edge]]
    end
    @assert QEM.left_face[cur_edge] == f1 
    @assert QEM.left_face[QEM.syms[cur_edge]] == f2
    return cur_edge 
end

@inline function UEdge(QEM::QuadEdgeMesh, e::QuadEdgeID) 
    UndirectedEdge(QEM.srcs[e], QEM.srcs[QEM.syms[e]])
end

function construct_pfaffian_orientation(
    edges::Set{UndirectedEdge}, nodes::Set{GlobalID}, 
    QEM::QuadEdgeMesh, rng, 
    edge_weights::Function 
)
    orientations = Dict{QuadEdgeID, Bool}()
    remaining_edges = edges 

    re = collect(remaining_edges)
    rev_vmap = collect(nodes) # maps 1...n to nodes 
    vmap = reverse_vmap(rev_vmap) #maps nodes to 1...n 
    if mod(length(vmap), 2) == 1 
        @assert false "attempting to make a pfaffian orientation on a graph with |V| odd" 
    end
    srcs, dsts = [vmap[e.id1] for e in re], [vmap[e.id2] for e in re]
    g = SimpleWeightedGraph(srcs, dsts, ones(length(re)))

    t = prim_mst(g)

    # randomly orient each edge in the tree 
    for edge in t 
        qe = QEM.edge_lookup[rev_vmap[edge.src], rev_vmap[edge.dst]]
        @assert typeof(qe) >: QuadEdgeID "no multiedges allowed when doing pfaffian"
        orientations[qe] = rand(rng, [true, false])
        orientations[QEM.syms[qe]] = !orientations[qe]
        delete!(remaining_edges, UndirectedEdge(rev_vmap[edge.src], rev_vmap[edge.dst]))
    end

    # construct the dual graph 
    remaining_edges_dual = Dict{QuadFaceID, Set{QuadEdgeID}}() 
    for edge in remaining_edges 
        QE = QEM.edge_lookup[(edge.id1, edge.id2)] 
        left_face = QEM.left_face[QE]
        right_face = QEM.left_face[QEM.syms[QE]]
        if haskey(remaining_edges_dual, left_face)
            push!(remaining_edges_dual[left_face], QE)
        else
            remaining_edges_dual[left_face] = Set([QE])
        end
        if haskey(remaining_edges_dual, right_face)
            push!(remaining_edges_dual[right_face], QEM.syms[QE])
        else
            remaining_edges_dual[right_face] = Set([QEM.syms[QE]])
        end
    end
    
    leaves_filtered = filter(((k,v),) -> length(v) == 1, remaining_edges_dual)
    leaves::Vector{Tuple{QuadFaceID, QuadEdgeID}} = [(k,first(v)) for (k,v) in leaves_filtered]

    i = 1

    #now, go through the leaves of the dual and pick the orientations of their edges 
    while length(remaining_edges_dual) > 1
        (leaf, og_edge) = leaves[i]
        nbr = QEM.left_face[QEM.syms[og_edge]]
        i += 1 
        nbr = first(nbr)
        @assert length(remaining_edges_dual[leaf]) == 1
        @assert QEM.left_face[og_edge] == leaf

        # count the number of clockwise_oriented edges around this face ("leaf")
        num_clockwise_around_face = 0 
        cur_edge = og_edge
        j = 0
        while true #loop around this face CCW 
            j += 1
            cur_edge = QEM.onext[QEM.syms[cur_edge]]
            @assert QEM.left_face[cur_edge] == leaf 
            if (cur_edge == og_edge) break end

            if !haskey(orientations, cur_edge)
                @show QEM.srcs[cur_edge] 
                @show leaf 
                @show QEM.srcs[og_edge]
                @show QEM.srcs[QEM.syms[og_edge]]
                @show haskey(orientations, og_edge)
            end 

            if orientations[cur_edge]
                num_clockwise_around_face += 1
            end 
        end

        # set the orientation of this edge 
        # want an odd number of clockwise total around each face 
        if mod(num_clockwise_around_face, 2) == 0 
            # make this one clockwise. 
            orientations[og_edge] = true 
            orientations[QEM.syms[og_edge]] = false
        else 
            orientations[og_edge] = false
            orientations[QEM.syms[og_edge]] = true 
        end 
        # delete this face ("leaf") from the dual
        delete!(remaining_edges_dual, leaf)
        delete!(remaining_edges_dual[nbr], QEM.syms[og_edge])
        if length(remaining_edges_dual[nbr]) == 1 
            #leaves[nbr] = remaining_edges_dual[nbr]
            push!(leaves, (nbr, first(remaining_edges_dual[nbr])))
        end
    end

    @assert length(leaves) == i

    rev_vmap = collect(nodes)#finest_level_nodes(G))
    vmap = reverse_vmap(rev_vmap)
    nv = length(rev_vmap)

    #now, make the tutte matrix using `orientations`
    for (_,face_edge) in QEM.faces 
        cur_edge = face_edge
        numcw = 0
        while true

            if !haskey(orientations, cur_edge) 
                @show face_edge 
                @show cur_edge 
                @show QEM.onext[QEM.syms[cur_edge]]
                @show QEM.onext[QEM.syms[face_edge]]
            end 

            if orientations[cur_edge]
                numcw += 1 
            end
            cur_edge = QEM.onext[QEM.syms[cur_edge]]
            if cur_edge == face_edge break end 
        end
        @assert mod(numcw, 2) == 1
    end

    tutte = zeros(nv, nv) #.+ (-1)
    for (k,v) in orientations
        if v == true 
            this_src = vmap[QEM.srcs[k]]
            this_dst = vmap[QEM.srcs[QEM.syms[k]]]
            x = edge_weights(QEM.srcs[k], QEM.srcs[QEM.syms[k]])
            tutte[this_src, this_dst] = x
            tutte[this_dst, this_src] = -x 
        end
    end
    return vmap, tutte
end

function log_count_perfect_matchings(nodes, edges, onbrs, rng)
    QEM = QuadEdgeMesh(nodes, onbrs)
    remove_multiedges!(QEM)
    vmap, tutte = construct_pfaffian_orientation(edges, nodes, QEM, rng, (_,_) -> 1)
    return 1/2 * logdet(tutte)
end

using InvertedIndices

function sample_perfect_matching(
    nodes::Set{GlobalID}, edges::Set{UndirectedEdge{GlobalID}}, 
    QEM::QuadEdgeMesh, 
    rng, 
    energy::Function # GlobalID * GlobalID -> Real 
)
    tutte = construct_pfaffian_orientation(nodes, edges, QEM, rng, energy)
    matching = Set{UndirectedEdge{GlobalID}}()

    G_nodes = finest_level_nodes(G)
    n = length(G_nodes)

    i = 0
    ms = []

    prev_num_matchings = 1/2 * logdet(tutte) 

    while length(matching) < n/2
        i += 1 
        if mod(i, 100) == 0
            push!(ms, length(matching))
            println(length(matching))
            return ms 
        end
        #if i == 100 
        #    return ms
        #end
        edge = findfirst(x->x != 0, tutte)
        @assert edge !== nothing 
        n1, n2 = Tuple(edge) 
        matchings_total = prev_num_matchings
        tutte[n1,n2] = 0 
        tutte[n2,n1] = 0
        matchings_without_edge = 1/2 * logdet(tutte) 

        p_edge_excluded = exp(matchings_without_edge - matchings_total)
        #take this edge with prob 1 - (matchings_without / matchings_total)
        #@show p_edge_excluded
        if rand(rng) <= 1 - p_edge_excluded
            push!(matching, UndirectedEdge(n1,n2)) 
            tutte = @view tutte[InvertedIndices.Not(n1, n2), InvertedIndices.Not(n1, n2)]
            # or in other words, tutte := tutte + tutte[n2, n1] * (e_n1 e_n2^T - e_n2 e_n1^T)
            prev_num_matchings = 1/2 * logdet(tutte)      
        else 
            prev_num_matchings = matchings_without_edge
            # we already zeroed that element of the tutte matrix. continue 
        end 
    end
    return ms, matching 
end

function sample_perfect_matching2(
    G, 
    QEM::QuadEdgeMesh, 
    rng, 
    energy::Function # GlobalID * GlobalID -> Real 
)
    tutte = construct_pfaffian_orientation(G, QEM, rng, energy)
    logdet_tutte = logdet(tutte)
    tutte_inv = inv(tutte)

    matching = Set{UndirectedEdge{GlobalID}}()
    n = length(finest_level_nodes(G))
    i = 0
    ms = []

    prev_num_matchings = 1/2 * logdet_tutte

    while length(matching) < n/2
        i += 1 
        if mod(i, 100) == 0
            push!(ms, length(matching))
            println(length(matching))
            return ms 
        end
        #if i == 100 
        #    return ms
        #end
        edge = findfirst(x->x != 0, tutte)
        @assert edge !== nothing 
        n1, n2 = Tuple(edge) 
        matchings_total = prev_num_matchings
        # really we don't need the 1/2 factor since it cancels but for clarity leave it for now
        
        tutte[n1,n2] = 0 
        tutte[n2,n1] = 0
        # or in other words, tutte := tutte + tutte[n2, n1] * (e_n1 e_n2^T - e_n2 e_n1^T)
        #logdet_tutte += log(1 + )

        matchings_without_edge = 1/2 * logdet(tutte) 

        p_edge_excluded = exp(matchings_without_edge - matchings_total)
        #take this edge with prob 1 - (matchings_without / matchings_total)
        #@show p_edge_excluded
        if rand(rng) <= 1 - p_edge_excluded
            push!(matching, UndirectedEdge(n1,n2)) 
            tutte = @view tutte[InvertedIndices.Not(n1, n2), InvertedIndices.Not(n1, n2)]
            prev_num_matchings = 1/2 * logdet(tutte)      
        else 
            prev_num_matchings = matchings_without_edge
            # we already zeroed that element of the tutte matrix. continue 
        end 
    end
    return ms, matching 
end

function tutte_remove_edge!(
    tutte::Matrix{Float64}, 
    tutte_inv::Matrix{Float64}, 
    logdet_tutte::Float64, i::Int, j::Int; 
    label=0
)    
    w = tutte[i,j]
    if w == 0 
        @assert tutte[j,i] == 0
        return logdet_tutte
    end

    logdet = add_symmetric_to_tutte!(tutte, tutte_inv, logdet_tutte, i, j, -w; label=label)
    #symmetrize_0s!(tutte_inv)
    return logdet
end

function tutte_remove_edge1!(
    tutte::Matrix{Float64}, 
    tutte_inv::Matrix{Float64}, 
    logdet_tutte::Float64, i::Int, j::Int; 
    label=0
)    
    w = tutte[i,j]
    if w == 0 
        @assert tutte[j,i] == 0
        return logdet_tutte
    end

    logdet1 = add_pointwise_to_tutte!(tutte, tutte_inv, logdet_tutte, i, j, -w; label=label)
    if logdet1 === nothing 
        return nothing 
    end
    logdet2 = add_pointwise_to_tutte!(tutte, tutte_inv, logdet1, j, i, w; label=label)
    if logdet2 === nothing 
        # need to put the first change back...
        logdet1 = add_pointwise_to_tutte!(tutte, tutte_inv, logdet_tutte, i, j, w; label=label)
        if logdet1 === nothing 
            @assert false "uh oh..." 
        end
        return nothing 
    end

    return logdet2 
end

function max_relative_deviation_from_skew_symmetric(tutte::Matrix{Float64}) 
    n,m = size(tutte) 
    max_d = 0 
    max_idx = (1,1)
    for i=1:n 
        for j=i:m
            d = relative_diff(tutte[i,j], tutte[j,i])
            if d > max_d 
                max_d = d 
                max_idx = (i,j)
            end
        end
    end
    return max_d, max_idx
end

function evaluate_inverseness(A, A_inv) 
    n,_ = size(A)
    return maximum(abs.(A * A_inv .- Matrix{Float64}(I, n,n)))
end

# updates tutte[i,j] += w and tutte[j,i] -= w. updates tutte_inv in place 
# and returns new logdet_tutte. 
function add_symmetric_to_tutte!(
    tutte::Matrix{Float64}, 
    tutte_inv::Matrix{Float64}, 
    logdet_tutte::Float64, i::Int, j::Int, w::Float64; 
    label="remove"
)
    #=n = size(tutte)[1]
    U = zeros(n,2); U[j, 1] = 1; U[i, 2] = 1
    V = zeros(2,n); V[1, i] = 1; V[2, j] = 1
    C = [-w 0; 0 w]=#

    #@assert tutte_inv[i,i] == 0 
    #@assert tutte_inv[j,j] == 0
    #@assert tutte_inv[j,i] == -tutte_inv[i,j]

    A_inv_U = tutte_inv[:, [j,i]]
    V_A_inv = tutte_inv[[i,j], :]
    C_inv_plus_V_A_inv_U::Matrix{Float64} = 
    [(tutte_inv[i,j] - 1/w) (tutte_inv[i,i]); (tutte_inv[j,j]) (tutte_inv[j,i] + 1/w)]

    det_factor_diff = det(C_inv_plus_V_A_inv_U) * (-w * w)
    det_factor_diff5 = (w * tutte_inv[i,j] - 1)^2

    if det_factor_diff < 1e-14
        #@show "is approx 0"
        return nothing 
    end

    inside_mat = inv(C_inv_plus_V_A_inv_U) * V_A_inv

    BLAS.gemm!('N', 'N', -1.0, A_inv_U, inside_mat, 1.0, tutte_inv)
    tutte[i,j] += w 
    tutte[j,i] -= w 

    #symmetrize_0s!(tutte_inv)
    return logdet_tutte + log(det_factor_diff5)
end

function add_pointwise_to_tutte!(tutte::Matrix{Float64}, 
    tutte_inv::Matrix{Float64}, 
    logdet_tutte::Float64, i::Int, j::Int, w::Float64; 
    label=""
)
    # u = e_i * w, v = e_j 
    if isapprox(w * tutte_inv[j,i], -1)
        #println("this decreases the determinant by a massive amount. add this edge to matching.")
        return nothing
    elseif w * tutte_inv[j,i] < -1 
        println("weird and negative, idk; label=", label)
        @show tutte_inv[j,i] 
        @show tutte_inv[i,j]
        @show w * tutte_inv[j,i]
        @show w * inv(tutte)[j,i]

        @show tutte_inv[j,i] - inv(tutte)[j,i]
        
        @show tutte[j,i], tutte[i,j]
        @assert false 
    else
        logdet_tutte += log(1 + w * tutte_inv[j,i])
        sc_fact = (-w/(1 + w * tutte_inv[j,i]))
        tutte[i,j] += w 
        # "GEneralized Matrix Multiplication"
        BLAS.gemm!('N', 'T', sc_fact, tutte_inv[:,i], tutte_inv[j,:], 1.0, tutte_inv)
        
        return logdet_tutte 
    end 
end

function relative_diff(a, b) 
    if abs(a) < abs(b) 
        return relative_diff(b,a) 
    end
    # a >= b 
    if isapprox(b, 0) 
        if isapprox(a, 0) return 0 
        else 
            if b == 0 
                return a 
            else 
                return (a-b)/b
            end
        end
    end
    return (a-b)/b 
end

# this one is actually slower.... 
function symmetrize_0s2!(A) 
    signs = A .< 0 
    A .+= -transpose(A) ./ 2
    A[A .< 0] .= 0 
    A[signs] .= -A[signs]
end

function symmetrize_0s!(A)
    n,m = size(A) 
    for i=1:n
        A[i,i] = 0 
        for j=i+1:m
            if A[i,j] == 0 || A[j,i] == 0 || (A[i,j] > 0 && A[j,i] > 0) || (A[i,j] < 0 && A[j,i] < 0)
                A[i,j] = 0 
                A[j,i] = 0 
            else 
                #they're different signs, as they should be 
                geomean = sqrt(-A[i,j] * A[j,i])
                if A[i,j] < 0 
                    A[i,j] = -geomean 
                    A[j,i] = geomean 
                else 
                    A[i,j] = geomean 
                    A[j,i] = -geomean 
                end
            end
        end
    end
end


# notable stuff for numerical stability: 
# - symmetrize_0s basically every time; knows about skew-symmetric 
# - skew symmetric edge deletion from adj matrix (fancy version of matrix det lemma)
#   (the woodbury identity)
# - explicitly casing out when the prob is very small, to avoid numerical stability Big Badness
# - check existence of PMs on unweighted ver of graph 

function fast_sample_pm(
    nodes::Set{GlobalID}, edges::Dict{GlobalID, Set{GlobalID}}, 
    QEM::QuadEdgeMesh, 
    rng, edge_weight::Function # GlobalID * GlobalID -> Real 
)

    edge_set = Set{UndirectedEdge}([UndirectedEdge(x,y) for (x, v) in edges for y in v])
    vmap, tutte::Matrix{Float64} = construct_pfaffian_orientation(edge_set, nodes, QEM, rng, edge_weight)

    return fast_sample_pm(nodes, edges, vmap, tutte, rng)
end

function fast_sample_pm(
    nodes::Set{GlobalID}, edges::Dict{GlobalID, Set{GlobalID}}, 
    vmap, tutte,
    rng # GlobalID * GlobalID -> Real 
)

    @assert all([length(v) !== 0 for (k,v) in edges])

    verbose = false #length(nodes) == 224

    matching = Set{UndirectedEdge{GlobalID}}()

    unweighted_tutte = copy(tutte) 
    for i=1:size(tutte)[1]
        for j=1:size(tutte)[1]
            unweighted_tutte[i,j] = tutte[i,j] == 0 ? 0 : (tutte[i,j] < 0 ? -1 : 1)
        end
    end
    try 
        log_unweighted_PMs = logdet(unweighted_tutte) 
        if log_unweighted_PMs < 0 
            println("det(unit. tutte) < 1; likely no PMs")
            return [], 0 
        end
        _ = inv(unweighted_tutte)
    catch 
        println("unit. tutte singular; likely no PMs")
        return [], 0
    end

    remaining_edges = deepcopy(edges)

    tutte_inv = []
    logdet_tutte = -1
    try 
        tutte_inv = inv(tutte) 
        logdet_tutte = logdet(tutte) 
    catch 
        #@assert false "unitary tutte matrix non-singular but weighted matrix singular??"
        println("weighted tutte matrix singular but unweighted non-singular; you likely have edge weights that are horrible which make the tutte matrix basically singular")
        return [], 0
    end

    p_incr = [] 
    log_p = 0

    i = 0
    num_nodes = length(nodes)
    #@show num_nodes 
    while length(matching) < num_nodes / 2 
        @assert all([length(v) !== 0 for (k,v) in remaining_edges])
        symmetrize_0s!(tutte_inv)
        i += 1

        if mod(i, 500) == 0 
            @show i, length(matching)
            tutte_inv = inv(tutte)
            symmetrize_0s!(tutte_inv)  
            logdet_tutte = logdet(tutte) 
        end

        e1, r = first(remaining_edges)
        e2 = first(r)
        delete!(remaining_edges[e1], e2)
        delete!(remaining_edges[e2], e1)

        if length(remaining_edges[e1]) == 0 || length(remaining_edges[e2]) == 0 
            if verbose 
                println("using by elim ", e1, " ", e2, " ", i)
            end 
            push!(matching, UndirectedEdge(e1, e2))
            for e1nbr in remaining_edges[e1]
                delete!(remaining_edges[e1nbr], e1)
            end 
            for e2nbr in remaining_edges[e2]
                delete!(remaining_edges[e2nbr], e2) 
            end
            delete!(remaining_edges, e1) 
            delete!(remaining_edges, e2)
            9
            @assert all([length(v) !== 0 for (k,v) in remaining_edges])
            continue 
        end

        n1, n2 = vmap[e1], vmap[e2]
        this_edge_weight = tutte[n1, n2]

        logdet_without_edge = tutte_remove_edge!(tutte, tutte_inv, logdet_tutte, n1, n2)

        if logdet_without_edge === nothing 
            # accept this edge 
            if verbose 
                println("using by necessity ", i)
            end 
            push!(matching, UndirectedEdge(e1, e2))
            for nbr in remaining_edges[e1]
                @assert nbr !== e2 
                delete!(remaining_edges[nbr], e1)
                logdet_tutte = tutte_remove_edge!(tutte, tutte_inv, logdet_tutte, vmap[e1], vmap[nbr];
                label="rem all besides necessary")
                if logdet_tutte === nothing
                    @show ":("
                    @assert false 
                end
            end 
            for nbr in remaining_edges[e2]
                delete!(remaining_edges[nbr], e2)
                logdet_tutte = tutte_remove_edge!(tutte, tutte_inv, logdet_tutte, vmap[e2], vmap[nbr];
                label=3)
                @assert logdet_tutte !== nothing
            end 
            delete!(remaining_edges, e1)
            delete!(remaining_edges, e2)
            @assert all([length(v) !== 0 for (k,v) in remaining_edges])
            continue 
        end

        p_excluded = exp(1/2 * logdet_without_edge - 1/2 * logdet_tutte)

        if !((p_excluded <= 1.0 && p_excluded >= 0.0) || isapprox(p_excluded, 1.0) || isapprox(p_excluded, 0.0))
            if isapprox(logdet_tutte, logdet_without_edge)
                # ok, whatever, rounding error i guess, and made worse by the exp(.) part 
            else 
                if p_excluded > 1 
                    println("p_excluded = ", p_excluded, " at i=", i, "; recomputing tutte inv and just dealing with it")
                    tutte_inv = inv(tutte) 
                else 
                    @show p_excluded, i
                    @assert false "p not between 0 and 1: ", p_excluded 
                end
            end
        end

        push!(p_incr, p_excluded)

        if rand(rng) < p_excluded 
            log_p += log(p_excluded)
            if verbose
                println("not using ", p_excluded, " ", e1, " ", e2, " ", i)
            end 
            @assert all([length(v) !== 0 for (k,v) in remaining_edges])
            # don't include this edge in the matching; delete it (we already did) and move on 
            logdet_tutte = logdet_without_edge 
        else 
            log_p += log(1-p_excluded)
            if verbose
                println("using ", p_excluded, " ", e1, " ", e2, " ", i) 
            end

            logdet_tutte = logdet_without_edge 

            # add this edge back 
            push!(matching, UndirectedEdge(e1, e2))
            logdet_tutte = add_symmetric_to_tutte!(tutte, tutte_inv, logdet_tutte, n1, n2, this_edge_weight > 0 ? 1.0 : -1.0; 
            label="add")
            @assert logdet_tutte !== nothing

            # delete all other incidences to n1 and n2. 
            for nbr in remaining_edges[e1]
                delete!(remaining_edges[nbr], e1)
                logdet_tutte = tutte_remove_edge!(tutte, tutte_inv, logdet_tutte, vmap[e1], vmap[nbr]; label="removenbrs")
                @assert logdet_tutte !== nothing
            end 
            for nbr in remaining_edges[e2]
                delete!(remaining_edges[nbr], e2)
                logdet_tutte = tutte_remove_edge!(tutte, tutte_inv, logdet_tutte, vmap[e2], vmap[nbr]; label="removenbrs")
                @assert logdet_tutte !== nothing
            end 
            @assert all([length(v) !== 0 for (k,v) in remaining_edges])

            delete!(remaining_edges, e1)
            delete!(remaining_edges, e2)
        end

        #if verbose 
        #    @assert isapprox(logdet_tutte, logdet(tutte))
        #end

    end
    @assert length(matching) == num_nodes / 2 
    return matching, log_p 
end

# edges in adjacency list form
function filter_edges!(edges, nodes) 
    for (k,v) in edges 
        if !(k in nodes)
            delete!(edges, k) 
        else
            edges[k] = filter(x -> x in nodes, v)
        end
    end
end

function even_ccs(nodes, edges) 
    conn_comps = Set{Set{GlobalID}}() 

    # basically do a DFS for each connected component 
    for k in nodes 
        this_cc = Set{GlobalID}()
        if any([k in x for x in conn_comps]) 
            continue 
        end
        frontier = Set{GlobalID}(k)
        while length(frontier) > 0 
            node = first(frontier)
            delete!(frontier, node)
            push!(this_cc, node)
            for nbr in edges[node]
                if nbr in this_cc || nbr in frontier continue end 
                push!(frontier, nbr)
            end
        end
        push!(conn_comps, this_cc)
    end
    #@show [length(x) for x in conn_comps]
    #@show count(x -> length(x) > 1, conn_comps)
    for cc in conn_comps
        if mod(length(cc), 2) == 1 
            this_cc_aps = articulation_points(cc, edges)
            removable_node = first(filter(x -> !(x in this_cc_aps), cc))
            delete!(cc, removable_node)
            if length(cc) == 0
                delete!(conn_comps, cc)
            end
        end
    end
    return conn_comps 
end
