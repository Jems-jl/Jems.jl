module StellarOpacity

using StellarChem, StellarConstants

abstract type AbstractOpacity end

include("SimpleElectronScattering.jl")

end # module StellarOpacity
