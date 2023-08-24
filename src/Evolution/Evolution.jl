module Evolution

export StellarModel

using Jems.Constants, Jems.Chem, Jems.EOS, Jems.Opacity

include("Options.jl")
include("StellarModel.jl")
include("Evaluation.jl")
include("Equations.jl")
include("IO.jl")
include("InitialCondition.jl")
include("EvolutionLoop.jl")

end # module Evolution
