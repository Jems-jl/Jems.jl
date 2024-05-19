module Jems

# each of these correspond to one of Jems' modules
include("DualSupport/DualSupport.jl")
include("Constants/Constants.jl")
include("Chem/Chem.jl")
include("EOS/EOS.jl")
include("ReactionRates/ReactionRates.jl")
include("NuclearNetworks/NuclearNetworks.jl")
include("Opacity/Opacity.jl")
include("Turbulence/Turbulence.jl")
include("StellarModels/StellarModels.jl")
include("Plotting/Plotting.jl")
include("Evolution/Evolution.jl")
include("DualExtrapolation/DualExtrapolation.jl") #new

end
