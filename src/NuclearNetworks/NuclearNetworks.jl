module NuclearNetworks

using ..ReactionRates, ..EOS

export NuclearNetwork, set_rates_for_network!, merge_nuclear_networks

abstract type AbstractNuclearNetwork end

@kwdef struct NuclearNetwork{TT} <: AbstractNuclearNetwork
    nspecies::Int  # Just the number of species in the network
    species_names::Vector{Symbol}  # just the species names
    reactions::TT
    xa_index::Dict{Symbol,Int}
    species_reactions_in::Vector{Vector{Tuple{Int,Int}}}
    species_reactions_out::Vector{Vector{Tuple{Int,Int}}}
end

function NuclearNetwork(species_names:: Vector{Symbol}, reactions::Vector{<: Any})
    nspecies = length(species_names)

    xa_index = Dict{Symbol, Int}()
    for (i,name) in enumerate(species_names)
        xa_index[name] = i
    end

    species_reactions_in = [Vector{Tuple{Int,Int}}() for _ in 1:nspecies]
    species_reactions_out = [Vector{Tuple{Int,Int}}() for _ in 1:nspecies]

    for (i, reaction) in enumerate(reactions)
        for j in eachindex(reaction.iso_in)
            iso_name = reaction.iso_in[j]
            if !haskey(xa_index, iso_name)
                error("Reaction uses $iso_name, but it is missing from species list.")
            end
            new_idx = xa_index[iso_name]
            stoich = reaction.num_iso_in[j]
            push!(species_reactions_in[new_idx], (i, stoich))
        end
        for j in eachindex(reaction.iso_out)
            iso_name = reaction.iso_out[j]
            if !haskey(xa_index, iso_name)
                error("Reaction uses $iso_name, but it is missing from species list.")
            end
            new_idx = xa_index[iso_name]
            stoich = reaction.num_iso_out[j]
            push!(species_reactions_out[new_idx], (i, stoich))
        end
    end

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

function NuclearNetwork(species_names::Vector{Symbol}, reaction_names::Vector{Tuple{Symbol, Symbol}})
    reactions = []
    

    for i in eachindex(reaction_names)
        reaction_name = reaction_names[i]
        reaction_obj = ReactionRates.reaction_list[reaction_name[1]][reaction_name[2]]
        push!(reactions, reaction_obj)
    end
    return NuclearNetwork(species_names, reactions)
end

function merge_nuclear_networks(nets::Vector{<:NuclearNetwork}) # merge_nuclear_networks([net1, net2, net3])
    all_species_list = Symbol[]
    all_reactions = []
    for net in nets
        append!(all_species_list, net.species_names)
        append!(all_reactions, net.reactions)
    end
    unique_species_names = unique(all_species_list)
    unique_reactions = unique(all_reactions)

    return NuclearNetwork(unique_species_names, unique_reactions)
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