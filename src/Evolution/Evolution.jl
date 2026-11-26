module Evolution

using Jems.Constants, Jems.Chem, Jems.EOS, Jems.Opacity, 
      Jems.NuclearNetworks, Jems.StellarModels, Jems.Plotting

include("Evaluation.jl")
include("Equations.jl")
include("LinearSolver.jl")
include("EvolutionInterface.jl")
include("NonlinearSolvers.jl")
include("StellarEvolution.jl")
include("OneZoneEvolution.jl")

end  # module Evolution
