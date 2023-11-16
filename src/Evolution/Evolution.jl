module Evolution

export StellarModel

using Jems.Constants, Jems.Chem, Jems.EOS, Jems.Opacity, 
      Jems.NuclearNetworks, Jems.StellarModels

include("Options.jl")
include("Plotter.jl")
include("StellarModel.jl")
include("Plotting.jl")
include("Evaluation.jl")
include("Equations.jl")
include("LinearSolver.jl")
include("EvolutionLoop.jl")

end  # module Evolution
