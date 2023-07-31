module StellarEvolution

export StellarModel

using StellarConstants
using StellarChem
using StellarEOS
using StellarOpacity

include("StellarModel.jl")
include("Solver.jl")
include("Equations.jl")

end # module StellarModel
