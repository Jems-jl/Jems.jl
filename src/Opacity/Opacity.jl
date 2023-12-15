module Opacity

using ..Chem, ..Constants

export AbstractOpacity

"""
    abstract type AbstractOpacity

Abstract supertype from which all defined opacity laws must derive, _ie_:

`struct MyOpacity <: AbstractOpacity`
"""
abstract type AbstractOpacity end

include("SimpleElectronScattering.jl")

end  # module Opacity
