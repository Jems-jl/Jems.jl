module Convection

using ..Constants

export AbstractConvection

abstract type AbstractConvection end

include("SimpleAdiabatic.jl")

end