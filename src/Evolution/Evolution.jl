module Evolution

export StellarModel

using Jems.Constants, Jems.Chem, Jems.EOS, Jems.Opacity, 
      Jems.NuclearNetworks, Jems.StellarModels

using LinearSolve

include("Evaluation.jl")
include("Equations.jl")
include("EvolutionLoop.jl")

end # module Evolution
