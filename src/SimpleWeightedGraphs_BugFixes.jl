""""""
function ne_bf(g::SimpleWeightedGraph)::Int
    return length(filter(e -> has_edge(g, e.src, e.dst), collect(edges(g))))
end

""""""
function connected_components_bf!(
    label::AbstractVector,
    g::SimpleWeightedGraph{T}
) where T
    nvg = nv(g)

    for u in vertices(g)
        label[u] != zero(T) && continue
        label[u] = u
        Q = Vector{T}()
        push!(Q, u)
        while !isempty(Q)
            src = popfirst!(Q)
            for vertex in all_neighbors(g, src)
                if g.weights[src, vertex] == 0
                    continue
                end
                if label[vertex] == zero(T)
                    push!(Q, vertex)
                    label[vertex] = u
                end
            end
        end
    end
    return label
end

""""""
function connected_components_bf(g::SimpleWeightedGraph{T}) where T
    label = zeros(T, nv(g))
    connected_components_bf!(label, g)
    c, d = Graphs.components(label)
    return c
end

""""""
function is_connected_bf(g::SimpleWeightedGraph)
    return ne_bf(g) + 1 >= nv(g) && length(connected_components_bf(g)) == 1
end