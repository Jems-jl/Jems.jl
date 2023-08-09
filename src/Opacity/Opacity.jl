module Opacity

using ..Chem, ..Constants

export AbstractOpacity

abstract type AbstractOpacity end

include("SimpleElectronScattering.jl")

end # module Opacity
