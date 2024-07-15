module ReactionRates

using ..Constants, ..EOS, ..Chem

export ToyReactionRate, KippReactionRate, JinaReactionRate, update_rate_cache!, RateCache

abstract type AbstractReactionRate end

reaction_list::Dict{Symbol, Dict{Symbol,AbstractReactionRate}} = Dict()

include("RateCache.jl")
include("JinaRates.jl")
include("KippRates.jl")
include("ToyRates.jl")

end