struct LevelSetToFractions{M, T, J, C}
    mesh::M
    topology::T
    jac::J
    colors::C
end
function LevelSetToFractions(mesh::Mesh)
    topo = topology(mesh)
    sparsity = get_sparsity_pattern(topo)
    jac = float.(sparsity)
    colors = matrix_colors(jac)
    return LevelSetToFractions(mesh, topo, jac, colors)
end
function Base.show(io::IO, ls::LevelSetToFractions)
    nelems = nelements(ls.mesh)
    nnodes = nvertices(ls.mesh)
    println(io, "Level set function: node values -> element volume fractions")
    println(io, "  #nodes = $nnodes \n  #elems = $nelems")
end
  
function get_sparsity_pattern(topo)
    I = Int[]
    J = Int[]
    V = Bool[]
    for (ei, e) in enumerate(elements(topo))
        for vi in indices(e)
            push!(I, ei)
            push!(J, vi)
            push!(V, true)
        end
    end
    return sparse(I, J, V, nelements(topo), nvertices(topo))
end

function (ls::LevelSetToFractions)(fracs::AbstractVector, C::AbstractVector)
    return volume_fractions!(fracs, ls.topology, ls.mesh, C)
end
function (ls::LevelSetToFractions)(C::AbstractVector)
    fracs = similar(C, nelements(ls.topology))
    return volume_fractions!(fracs, ls.topology, ls.mesh, C)
end
function ChainRulesCore.rrule(ls::LevelSetToFractions, C::AbstractVector)
    forwarddiff_color_jacobian!(ls.jac, ls, C, colorvec = ls.colors)
    return ls(C), Δ -> begin
        (NoTangent(), ls.jac' * Δ)
    end
end
