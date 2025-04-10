"""
    struct ToyReactionRate{TT<:Real}<:ReactionRates.AbstractReactionRate


- `name`: `Symbol` giving the name of the reaction.
- `iso_in`: vector that contains all reactants given as symbols (e.g. `[:H1, :H2]`)
- `num_iso_in`: number of each of the elements in `iso_in` that are used in the reaction, given as a vector of integers. For example if `iso_in` is `[:He4]` and `num_iso_in` is `[3]` it means the reaction uses three ":He4".
- `iso_out`: Same as `iso_in` but for the products of the reaction.
- `num_iso_out`: Same as `num_iso_in` but for the products of the reaction.
- `Qvalue`: Q-value of the reaction (in erg), simply given by the mass difference.
"""
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

"""
    function get_reaction_rate(reaction::ToyReactionRate, T::T1, ρ::T2, xa::AbstractVector{TT},
                                xa_index::Dict{Symbol,Int})::TT where {TT,T1,T2}

Input:
reaction: the reaction to evaluate for
- `T``: the temperature in K
- `ρ``: the density g cm^-3
- `xa``: element mass fractions
- `xa_index`: Dictionary containing the index of each element within `xa`

Output:
- ϵ_nuc / Qvalue, has units s^-1 g^-1
"""
function get_reaction_rate(reaction::ToyReactionRate, T::T1, ρ::T2, xa::AbstractVector{TT},
                           xa_index::Dict{Symbol,Int})::TT where {TT,T1,T2}
    if reaction.name == :toy_pp  # Taken from Onno Pols' lectures with arbitrary pre-factor
        ϵnuc = 0.1 * xa[xa_index[:H1]]^2 * ρ * (T / 1e6)^4
        return ϵnuc / reaction.Qvalue
    elseif reaction.name == :toy_cno  # Taken from Onno Pols' lectures with arbitrary pre-factor
        ϵnuc = 0.1 * xa[xa_index[:H1]]^2 * ρ * (T / 1e7)^18
        return ϵnuc / reaction.Qvalue
    else
        throw(ArgumentError("No method to compute rate for $(reaction.name)"))
    end
end
