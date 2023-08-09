module Opacity

using Jems.Chem, Jems.Constants

export AbstractOpacity

abstract type AbstractOpacity end

include("SimpleElectronScattering.jl")

end # module Opacity
