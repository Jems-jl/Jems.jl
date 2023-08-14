function simple_adiabatic(sm::StellarModel, k::Int,
                          varm1::Vector{<:TT}, var00::Vector{<:TT}, varp1::Vector{<:TT},
                          eosm1::Vector{<:TT}, eos00::Vector{<:TT}, eosp1::Vector{<:TT},
                          κm1::TT, κ00::TT, κp1::TT)::TT where {TT<:Real}
    return (sm.dm[k] * eos00[7] + sm.dm[k + 1] * eosp1[7]) / (sm.dm[k] + sm.dm[k + 1])
end