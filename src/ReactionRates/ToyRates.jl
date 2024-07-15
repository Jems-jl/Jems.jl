struct ToyReactionRate{TT<:Real} <: AbstractReactionRate
    name::Symbol
    iso_in::Vector{Symbol}
    num_iso_in::Vector{Int64}
    iso_out::Vector{Symbol}
    num_iso_out::Vector{Int64}
    Qvalue::TT
end

reaction_list[:toy_rates] = Dict(
    :toy_pp => ToyReactionRate(:toy_pp, [:H1], [4], [:He4], [1],
                               ((4 * Chem.isotope_list[:H1].mass - Chem.isotope_list[:He4].mass) * AMU * CLIGHT^2)),
    :toy_cno => ToyReactionRate(:toy_cno, [:H1], [4], [:He4], [1],
                                ((4 * Chem.isotope_list[:H1].mass - Chem.isotope_list[:He4].mass) * AMU * CLIGHT^2)),
)

function get_reaction_rate(reaction::ToyReactionRate, cache::RateCache{TT}, ρ::TT, xa::AbstractVector{TT},
                           xa_index::Dict{Symbol,Int})::TT where {TT}
    if reaction.name == :toy_pp  # Taken from Onno Pols' lectures with arbitrary pre-factor
        ϵnuc = 0.1 * xa[xa_index[:H1]]^2 * ρ * (cache.T9 / 1e3)^4
        return ϵnuc / reaction.Qvalue
    elseif reaction.name == :toy_cno  # Taken from Onno Pols' lectures with arbitrary pre-factor
        ϵnuc = 0.1 * xa[xa_index[:H1]]^2 * ρ * (cache.T9 / 1e2)^18
        return ϵnuc / reaction.Qvalue
    else
        throw(ArgumentError("No method to compute rate for $(reaction.name)"))
    end
end
