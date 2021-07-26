function volume_fractions(topo, mesh, C)
    map(1:nelements(mesh)) do i
        el = element(mesh, i)
        T = typeof(e)
        inds = Tuple(indices(element(topo, i)))
        vs = Tuple(vertices(el))
        Cs = C[inds]
        return positive_measure_ratio(T, inds, vs, Cs)
    end
end

function positive_measure_ratio(::Type{<:Triangle}, inds, vs, Cs)
    a1, a2 = _positive_measures(Triangle, inds, vs, Cs)
    return a1 / a2
end
function _positive_measures(::Type{<:Triangle}, inds, vs, Cs)
    nneg = 0
    npos = 0
    all_pos = true
    all_neg = true
    for ind in 1:3
        if Cs[ind] < 0
            all_pos = false
            nneg += 1
        end
        if Cs[ind] > 0
            all_neg = false
            npos += 1
        end
    end
    if all_neg
        return zero(eltype(Cs))
    elseif all_pos
        return one(eltype(Cs))
    end

    ind0 = -1
    for j in 1:3
        if nneg == 1 && C[inds[j]] < 0
            ind0 = j
            break
        elseif nneg == 2 && C[inds[j]] > 0
            ind0 = j
            break
        end
    end
    @assert ind0 != -1
    others = setdiff((1, 2, 3), ind0)

    δC1 = Cs[ind0] - Cs[others[1]]
    δx1 = vs[ind0] - vs[others[1]]
    ratio1 = Cs[ind0] / δC1 * norm(δx1)
    x1 = vs[ind0] - ratio1 * δx1

    δC2 = Cs[ind0] - Cs[others[2]]
    δx2 = vs[ind0] - vs[others[2]]
    ratio2 = Cs[ind0] / δC2 * norm(δx2)
    x2 = vs[ind0] - ratio2 * δx2

    a1 = Meshes.measure(Triangle(vs...)) |> abs
    a2 = Meshes.measure(Triangle(x1, x2, vs[ind0])) |> abs
    
    if nneg == 1
        return (a1 - a2), a1
    else
        return a2, a1
    end
end

function positive_measure_ratio(::Type{<:Quadrangle}, inds, vs, Cs)
    a1, a2 = _positive_measures(Quadrangle, inds, vs, Cs)
    return a1 / a2
end
function _positive_measures(::Type{<:Quadrangle}, inds, vs, Cs)
    inds1 = @SVector [1, 2, 3]
    a1, b1 = _positive_measures(Triangle, inds[inds1], vs[inds1], Cs[inds1])
    inds2 = @SVector [1, 3, 4]
    a2, b2 = _positive_measures(Triangle, inds[inds2], vs[inds2], Cs[inds2])
    return (a1 + a2), (b1 + b2)
end

function positive_measure_ratio(::Type{<:Tetrahedron}, inds, vs, Cs)
    a1, a2 = _positive_measures(Tetrahedron, inds, vs, Cs)
    return a1 / a2
end
function _positive_measures(::Type{<:Tetrahedron}, inds, vs, Cs)
    nneg = 0
    npos = 0
    all_pos = true
    all_neg = true
    for ind in 1:4
        if Cs[ind] < 0
            all_pos = false
            nneg += 1
        end
        if Cs[ind] > 0
            all_neg = false
            npos += 1
        end
    end
    if all_neg
        return zero(eltype(Cs))
    elseif all_pos
        return one(eltype(Cs))
    end

    neg_inds = Int[]
    for j in 1:4
        if C[inds[j]] < 0
            push!(neg_inds, j)
        end
    end
    if nneg == 1
        others = setdiff((1, 2, 3, 4), neg_inds)
        ind0 = neg_inds[1]

        δC1 = Cs[ind0] - Cs[others[1]]
        δx1 = vs[ind0] - vs[others[1]]
        ratio1 = Cs[ind0] / δC1 * norm(δx1)
        x1 = vs[ind0] - ratio1 * δx1
    
        δC2 = Cs[ind0] - Cs[others[2]]
        δx2 = vs[ind0] - vs[others[2]]
        ratio2 = Cs[ind0] / δC2 * norm(δx2)
        x2 = vs[ind0] - ratio2 * δx2

        δC3 = Cs[ind0] - Cs[others[3]]
        δx3 = vs[ind0] - vs[others[3]]
        ratio3 = Cs[ind0] / δC3 * norm(δx3)
        x3 = vs[ind0] - ratio3 * δx3

        a1 = Meshes.measure(Tetrahedron(vs...)) |> abs
        a2 = Meshes.measure(Tetrahedron(x1, x2, x3, vs[ind0])) |> abs

        return (a1 - a2), a1
    elseif length(neg_inds) == 2
        others = setdiff((1, 2, 3, 4), neg_inds)
        ind0, ind1 = neg_inds

        δC1 = Cs[ind0] - Cs[others[1]]
        δx1 = vs[ind0] - vs[others[1]]
        ratio1 = Cs[ind0] / δC1 * norm(δx1)
        x1 = vs[ind0] - ratio1 * δx1
    
        δC2 = Cs[ind0] - Cs[others[2]]
        δx2 = vs[ind0] - vs[others[2]]
        ratio2 = Cs[ind0] / δC2 * norm(δx2)
        x2 = vs[ind0] - ratio2 * δx2

        δC3 = Cs[ind1] - Cs[others[1]]
        δx3 = vs[ind1] - vs[others[1]]
        ratio3 = Cs[ind1] / δC3 * norm(δx3)
        x3 = vs[ind1] - ratio3 * δx3

        δC4 = Cs[ind1] - Cs[others[2]]
        δx4 = vs[ind1] - vs[others[2]]
        ratio4 = Cs[ind1] / δC4 * norm(δx4)
        x4 = vs[ind1] - ratio4 * δx4

        a1 = Meshes.measure(Tetrahedron(vs...)) |> abs
        a2 = Meshes.measure(Wedge(x1, x2, vs[ind0], x3, x4, vs[ind1])) |> abs

        return (a1 - a2), a1
    elseif length(neg_inds) == 3
        pos_inds = setdiff((1, 2, 3, 4), neg_inds)
        ind0 = pos_inds[1]
        others = neg_inds

        δC1 = Cs[ind0] - Cs[others[1]]
        δx1 = vs[ind0] - vs[others[1]]
        ratio1 = Cs[ind0] / δC1 * norm(δx1)
        x1 = vs[ind0] - ratio1 * δx1
    
        δC2 = Cs[ind0] - Cs[others[2]]
        δx2 = vs[ind0] - vs[others[2]]
        ratio2 = Cs[ind0] / δC2 * norm(δx2)
        x2 = vs[ind0] - ratio2 * δx2

        δC3 = Cs[ind0] - Cs[others[3]]
        δx3 = vs[ind0] - vs[others[3]]
        ratio3 = Cs[ind0] / δC3 * norm(δx3)
        x3 = vs[ind0] - ratio3 * δx3

        a1 = Meshes.measure(Tetrahedron(vs...)) |> abs
        a2 = Meshes.measure(Tetrahedron(x1, x2, x3, vs[ind0])) |> abs

        return a2, a1
    else
        throw("Unreachable.")
    end
end

function positive_measure_ratio(::Type{<:Hexahedron}, inds, vs, Cs)
    a1, a2 = _positive_measures(Hexahedron, inds, vs, Cs)
    return a1 / a2
end
function _positive_measures(::Type{<:Hexahedron}, inds, vs, Cs)
    _inds = SVector(1, 2, 4, 5)
    a1, b1 = _positive_measures(Tetrahedron, inds[_inds], vs[_inds], Cs[_inds])
    _inds = SVector(2, 7, 4, 5)
    a2, b2 = _positive_measures(Tetrahedron, inds[_inds], vs[_inds], Cs[_inds])
    _inds = SVector(4, 5, 7, 8)
    a3, b3 = _positive_measures(Tetrahedron, inds[_inds], vs[_inds], Cs[_inds])
    _inds = SVector(2, 6, 7, 5)
    a4, b4 = _positive_measures(Tetrahedron, inds[_inds], vs[_inds], Cs[_inds])
    _inds = SVector(4, 7, 2, 3)
    a5, b5 = _positive_measures(Tetrahedron, inds[_inds], vs[_inds], Cs[_inds])
    return (a1 + a2 + a3 + a4 + a5), (b1 + b2 + b3 + b4 + b5)
end
