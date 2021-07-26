"""
    Wedge(p1, p2, p3, p4, p5, p6)

A wedge with points `p1`, `p2`, `p3`, `p4`, `p5`, `p6`.
"""
struct Wedge{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Meshes.Polyhedron{Dim,T}
  vertices::V
end

Meshes.nvertices(::Type{<:Wedge}) = 6
Meshes.nvertices(p::Wedge) = nvertices(typeof(p))

Meshes.isconvex(::Type{<:Wedge}) = true
Meshes.issimplex(::Type{<:Wedge}) = false

function Meshes.measure(t::Wedge)
  # Triangle 1 - side (1, 2, 3)
  a1 = t.vertices[1].coords
  n1 = (t.vertices[2].coords - t.vertices[1].coords) × (t.vertices[3].coords - t.vertices[1].coords)

  # Quad 1 - side (1, 3, 6, 4)
  a2 = t.vertices[1].coords
  n2 = (t.vertices[3].coords - t.vertices[1].coords) × (t.vertices[6].coords - t.vertices[1].coords) + 
       (t.vertices[6].coords - t.vertices[1].coords) × (t.vertices[4].coords - t.vertices[1].coords)

  # Triangle 2 - side (4, 6, 5)
  a3 = t.vertices[4].coords
  n3 = (t.vertices[6].coords - t.vertices[4].coords) × (t.vertices[5].coords - t.vertices[4].coords)

  # Quad 2 - side (2, 5, 6, 3)
  a4 = t.vertices[2].coords
  n4 = (t.vertices[5].coords - t.vertices[2].coords) × (t.vertices[6].coords - t.vertices[2].coords) + 
       (t.vertices[6].coords - t.vertices[2].coords) × (t.vertices[3].coords - t.vertices[2].coords)

  # Quad 3 - base (1, 4, 5, 2)
  a5 = t.vertices[1].coords
  n5 = (t.vertices[4].coords - t.vertices[1].coords) × (t.vertices[5].coords - t.vertices[1].coords) + 
       (t.vertices[5].coords - t.vertices[1].coords) × (t.vertices[2].coords - t.vertices[1].coords)
  
  s = dot(a1, n1) + dot(a2, n2) + dot(a3, n3) + dot(a4, n4) + dot(a5, n5)
  return s / 6
end
