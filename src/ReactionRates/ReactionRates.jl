module ReactionRates

using ..Constants, ..EOS, ..Chem

export ToyReactionRate, KippReactionRate

abstract type AbstractReactionRate end

reaction_list::Dict{Symbol, Dict{Symbol,AbstractReactionRate}} = Dict()

include("JinaRates.jl")
include("KippRates.jl")
include("ToyRates.jl")

end