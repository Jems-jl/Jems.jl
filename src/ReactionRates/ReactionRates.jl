module ReactionRates

using ..Constants, ..EOS, ..Chem

export ToyReactionRate, KippReactionRate

"""
    abstract type AbstractReactionRate end

Reaction rates should be defined as a subtype of `AbstractReactionRate`. Any custom reaction rate needs to contain a set of fixed fields:

```julia
struct MyReactionRate{TT<:Real} <: ReactionRates.AbstractReactionRate
    name::Symbol
    # num_iso_in of isotopes iso_in are converted into num_iso_out of isotopes iso_out
    iso_in::Vector{Symbol}
    num_iso_in::Vector{Int}
    iso_out::Vector{Symbol}
    num_iso_out::Vector{Int}
    Qvalue::TT  # energy released per reaction of this type (i.e. the mass defect), in erg
end
```
arbitrary fields can be added here. For instance if you have a tabulated reaction rate you can add the table values or even include the data defining a specific interpolator for it.

For any custom rate the following function needs to be defined
```julia
function get_reaction_rate(reaction::MyReactionRate, T::T1, ρ::T2, xa::AbstractVector{TT}
    #insert your code here
end
```
the function computes the specific reaction rate based on the `name` field of the reaction or done generically based on the other fields of the struct. Units are in ``\\mathrm{s^{-1}\\;g^{-1}}``
"""
abstract type AbstractReactionRate end

"""
    reaction_list::Dict{Symbol, Dict{Symbol,AbstractReactionRate}} = Dict()

The `reaction_list` dictionary contains reaction rates that are available for construction of nuclear reaction networks. It is a dictionary of dictionaries, for instance accessing `reaction_list[:jina_rates]` will provide all available rates from JINA.

You can add a new dictionary of rates to `reaction_list` by doing
```julia
reaction_list[:kipp_rates] = Dict(
    :kipp_pp => MyReactionRate(:my_pp, [:H1], [2], [:H2], [1],
                                 ((2 * Chem.isotope_list[:H1].mass - Chem.isotope_list[:He2].mass) * AMU * CLIGHT^2)),
                                 # add whatever more rates
                                 )
```
this assumes a specific form for the constructor of `MyReactionRate` which only includes the mandatory fields, you can use any custom constructor for rate definitions which have more fields than the standard ones.
"""
reaction_list::Dict{Symbol, Dict{Symbol,AbstractReactionRate}} = Dict()

include("JinaRates.jl")
include("KippRates.jl")
include("ToyRates.jl")

end