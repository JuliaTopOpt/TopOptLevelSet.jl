module TopOptLevelSet

using TupleTools, Meshes, StaticArrays
const TOLS = TopOptLevelSet

export TOLS

include("wedge.jl")
include("fractions.jl")

end
