module ReactionRates

using ..Constants, ..EOS, ..Chem

export ToyReactionRate, KippReactionRate

abstract type AbstractReactionRate end

"""

    const reaction_list::Dict{Symbol, Dict{Symbol,AbstractReactionRate}} 

Dictionary that contains all reaction rates  currently implemented in Jems.
To access a rate, write, _e.g.:_

    reaction_list[:kipp_rates][:kipp_pp]

"""
const reaction_list::Dict{Symbol, Dict{Symbol,AbstractReactionRate}} = Dict()

include("KippRates.jl")
include("ToyRates.jl")

end