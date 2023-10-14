module StellarModels

export StellarModel, get_tmp

using Jems.Constants, Jems.Chem, Jems.EOS, Jems.Opacity, Jems.NuclearNetworks

include("Options.jl")
include("StellarModel.jl")
include("IO.jl")
include("InitialCondition.jl")
include("Remesher.jl")

end # module Evolution
