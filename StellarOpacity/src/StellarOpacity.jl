module StellarOpacity

using StellarChem, StellarConstants

export AbstractOpacity

abstract type AbstractOpacity end

include("SimpleElectronScattering.jl")

end # module StellarOpacity
