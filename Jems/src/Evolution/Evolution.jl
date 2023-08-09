module Evolution

export StellarModel

using Jems.Constants, Jems.Chem, Jems.EOS, Jems.Opacity

include("Options.jl")
include("StellarModel.jl")
include("Solver.jl")
include("Equations.jl")
include("IO.jl")
include("InitialCondition.jl")
include("EvolutionLoop.jl")

end # module Evolution
