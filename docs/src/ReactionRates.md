# ReactionRates

The [ReactionRates](@ref) module of Jems handles all things related to the nuclear reaction rates of mixtures.

## General functionality
The [`ReactionRates`](@ref) module has a simple structure, rates can be defined as subtypes of the abstract type [`Jems.ReactionRates.AbstractReactionRate`](@ref). A function needs to be defined to compute the rate which dispatches based on the type of the reaction rate. Multiple rates can also be added to the `Jems.reactions_list` dictionary in order to construct nuclear reaction networks with the `NuclearNetworks` module.

A simple example for this would be:

```@example
using Jems.ReactionRates
using Jems.Chem
using Jems.Constants

struct MyReactionRate{TT<:Real} <: ReactionRates.AbstractReactionRate
    name::Symbol
    iso_in::Vector{Symbol}
    num_iso_in::Vector{Int}
    iso_out::Vector{Symbol}
    num_iso_out::Vector{Int}
    Qvalue::TT
end

function ReactionRates.get_reaction_rate(reaction::MyReactionRate, T::T1, ρ::T2, xa::AbstractVector{TT}, xa_index::Dict{Symbol,Int})::TT where {TT,T1,T2}
    rate = 0
    if reaction.name == :my_h1_h1_to_d2
        h1_index = xa_index[:H1]
        h1 = xa[h1_index]
        rate = h1^2 # mock example, we just have h1^2 reactions per second per gram.
    end
    return rate
end

my_reaction = MyReactionRate(:my_h1_h1_to_d2, [:H1], [2], [:D2], [1], (2 * Chem.isotope_list[:H1].mass - Chem.isotope_list[:D2].mass) * AMU * CLIGHT^2)

# evaluate this rate
T = 1e9 # in K
ρ = 10 # in g cm^-3
xa = [1.0, 0.0] # mass fractions, considering only H1 and H2, with 100% H1
xa_index = Dict(:H1=>1, :H2=>2)
ReactionRates.get_reaction_rate(my_reaction, T, ρ, xa, xa_index)

# store in reactions_list, this is done by creating a new dictionary for my rates, in this case with a single entry
ReactionRates.reaction_list[:my_rates] = Dict(:my_h1_h1_to_d2 => my_reaction)
;
```
For more details check the documentation of [`Jems.ReactionRates.AbstractReactionRate`](@ref) and [`Jems.ReactionRates.reaction_list`](@ref).

```@docs
Jems.ReactionRates.AbstractReactionRate
Jems.ReactionRates.reaction_list
```

## Kippenhahn rates
The Kippenhahn reaction rates are based on the [Kippenhahn textbook](https://doi.org/10.1007/978-3-642-30304-3) on stellar structure and evolution. Rates are estimated simply by taking the values of $\epsilon_\mathrm{nuc}$ derived using the formulae in the textbook and dividing by the ``Q`` value of the reaction, taken to be simply the mass difference. Rates available are:
- `:kipp_pp`: compound rate for the full pp-chain, takes 4 H1 and makes one He4
- `:kipp_cno`: compund rate for the full cno-chain, takes 4 H1 and makes one He4
- `:kipp_3alphaCF88`: TODO
- `:kipp_3alphaA99`: TODO
- `:kipp_C12alpha`: TODO
- `:kipp_O16alpha`: TODO
- `:kipp_CC`: TODO
- `:kipp_OO`: TODO
The rates are stored in `reactions_list[:kipp_rates]`, so they can be accesed as, for example, `reactions_list[:kipp_rates][:kipp_pp]`.


```@docs
Jems.ReactionRates.KippReactionRate
Jems.ReactionRates.get_reaction_rate(reaction::Jems.ReactionRates.KippReactionRate, T::T1, ρ::T2, xa::AbstractVector{TT},xa_index::Dict{Symbol,Int}) where{TT,T1,T2}
```

## JINA rates
TODO: Jina rates have a bunch of extra functionality to help access independent reactions from their large library, need to describe these and possibly clean up the code there a bit.
```@docs
Jems.ReactionRates.JinaReactionRate
```

## Toy rates

Early rates used within the code based on the power law relationships given in the lecture notes on stellar evolution of Onno Pols. The lectures only provided relationships of the form ``\\epsilon_\\mathrm{nuc}\\propto \rho T^\\nu``, here we made an arbitrary choice for the pre-factor and determine the rate by dividing for the ``Q`` value. **YOU ARE STRONGLY ADVISED NOT TO USE THESE RATES**.

Rates available are:
- `:toy_pp`: Compound rate for the pp-chain. Uses $\epsilon_\mathrm{nuc}\propto\rho T^4$.
- `:toy_cno`: Compound rate for the CNO cycle. Uses $\epsilon_\mathrm{nuc}\propto\rho T^{18}$.
The rates are stored in `reactions_list[:toy_rates]`, so they can be accesed as, for example, `reactions_list[:kipp_rates][:toy_pp]`.

```@docs
Jems.ReactionRates.ToyReactionRate
Jems.ReactionRates.get_reaction_rate(reaction::Jems.ReactionRates.ToyReactionRate, T::T1, ρ::T2, xa::AbstractVector{TT},xa_index::Dict{Symbol,Int}) where{TT,T1,T2}
```