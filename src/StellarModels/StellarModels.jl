module StellarModels

export AbstractModel, AbstractModelProperties, OneZone, OneZoneProperties, StellarModel, StellarModelProperties

using Jems.Constants, Jems.Chem, Jems.EOS, Jems.Opacity, Jems.NuclearNetworks

abstract type AbstractModel end
abstract type AbstractModelProperties end

include("Options.jl")
include("PlotterInterface.jl")
include("StellarModelProperties.jl")
include("SolverData.jl")
include("EquationSupport.jl")
include("StellarModel.jl")
include("OneZone.jl")
include("IO.jl")
include("InitialCondition.jl")
include("Remesher.jl")

end  # module StellarModels
