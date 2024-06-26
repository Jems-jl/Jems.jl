module NuclearNetworks

using ..ReactionRates, ..EOS

export NuclearNetwork, set_rates_for_network!

abstract type AbstractNuclearNetwork end

@kwdef struct NuclearNetwork{TT} <: AbstractNuclearNetwork
    nspecies::Int  # Just the number of species in the network
    species_names::Vector{Symbol}  # just the species names
    reactions::TT
    xa_index::Dict{Symbol,Int}
    species_reactions_in::Vector{Vector{Tuple{Int,Int}}}
    species_reactions_out::Vector{Vector{Tuple{Int,Int}}}
end

function NuclearNetwork(species_names, reaction_names::Vector{Tuple{Symbol, Symbol}})
    nspecies = length(species_names)
    xa_index = Dict{Symbol, Int}()
    for i in eachindex(species_names)
        xa_index[species_names[i]] = i
    end

    reactions = []
    species_reactions_in::Vector{Vector{Tuple{Int,Int}}} = [[] for i in eachindex(species_names)]
    species_reactions_out::Vector{Vector{Tuple{Int,Int}}} = [[] for i in eachindex(species_names)]
    for i in eachindex(reaction_names)
        reaction_name = reaction_names[i]
        reaction = ReactionRates.reaction_list[reaction_name[1]][reaction_name[2]]
        push!(reactions, reaction)
        for j in eachindex(reaction.iso_in)
            if ! (reaction.iso_in[j] ∈ species_names)
                throw(ArgumentError("Reaction $(reaction.name) requires isotope $(reaction.iso_in[j]), but it is not part of species_names"))
            end
            push!(species_reactions_in[xa_index[reaction.iso_in[j]]], (i,reaction.num_iso_in[j]))
        end
        for j in eachindex(reaction.iso_out)
            if ! (reaction.iso_out[j] ∈ species_names)
                throw(ArgumentError("Reaction $(reaction.name) requires isotope $(reaction.iso_out[j]), but it is not part of species_names"))
            end
            push!(species_reactions_out[xa_index[reaction.iso_out[j]]], (i,reaction.num_iso_out[j]))
        end
    end

    # At this point reactions is a Vector{Any}. For type stability we want to turn it into a vector with type
    # Union{...}, where the Union contains all types of reactions

    reactions_typed::Vector{Union{typeof.(reactions)...}} = [reactions...]

    return NuclearNetwork(
        nspecies = nspecies,
        species_names = species_names,
        reactions = reactions_typed,
        xa_index = xa_index,
        species_reactions_in = species_reactions_in,
        species_reactions_out = species_reactions_out,
    )
end

function set_rates_for_network!(rates::AbstractArray{TT}, net::NuclearNetwork, T::T1, ρ::T2,
                                xa::AbstractArray{TT}) where {TT,T1,T2}
    if length(rates) != length(net.reactions)
        throw(ArgumentError("Length of `rates` and `net.reactions` must be equal"))
    end
    for i in eachindex(rates)
        rates[i] = ReactionRates.get_reaction_rate(net.reactions[i], T, ρ, xa, net.xa_index)
    end
end

end