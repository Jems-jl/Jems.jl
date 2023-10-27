module NuclearNetworks

using ..ReactionRates, ..EOS

export NuclearNetwork, set_rates_for_network!

@kwdef struct NuclearNetwork{TT}
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

    reactions::Vector{Union{ToyReactionRate{Float64}}} = []
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

    return NuclearNetwork(
        nspecies = nspecies,
        species_names = species_names,
        reactions = reactions,
        xa_index = xa_index,
        species_reactions_in = species_reactions_in,
        species_reactions_out = species_reactions_out,
    )
end

function set_rates_for_network!(rates::AbstractArray{TT}, net::NuclearNetwork, eos00::EOSResults{TT}, xa::AbstractArray{TT}) where{TT}
    if length(rates) != length(net.reactions)
        throw(ArgumentError("Length of `rates` and `net.reactions` must be equal"))
    end
    for i in eachindex(rates)
        rates[i] = ReactionRates.get_reaction_rate(net.reactions[i], eos00, xa, net.xa_index)
    end
    return
end

end