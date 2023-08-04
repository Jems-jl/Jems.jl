module StellarEvolution

export StellarModel

using StellarConstants
using StellarChem
using StellarEOS
using StellarOpacity

include("Options.jl")
include("StellarModel.jl")
include("Solver.jl")
include("Equations.jl")
include("IO.jl")
#include("EvolutionLoop.jl")

end # module StellarModel
