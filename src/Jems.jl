module Jems

# each of these correspond to one of Jems' modules
include("Constants/Constants.jl")
include("Chem/Chem.jl")
include("EOS/EOS.jl")
include("ReactionRates/ReactionRates.jl")
include("NuclearNetworks/NuclearNetworks.jl")
include("Opacity/Opacity.jl")
include("StellarModels/StellarModels.jl")
include("Evolution/Evolution.jl")

end
