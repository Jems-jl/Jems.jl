export get_opacity_resultsTρ

struct KramerPrescripted <: AbstractOpacity end

"""
    get_opacity_resultsTP(opacity::SimpleElectronScatteringOpacity, lnT::TT, lnP::TT, xa::Vector{<:TT},
                            species::Vector{Symbol})::TT where {TT<:Real}

Evaluates the opacity of the current mixture with mass fractions `xa`, species symbols `species` (both these should be
of length `nspecies`), the natural log of temperature and density `lnT`, `lnρ`, and the opacity law `opacity`.
"""
function get_opacity_resultsTρ(opacity::KramerPrescripted, lnT::TT, lnρ::TT,
                               xa::AbstractVector{<:TT}, species::Vector{Symbol})::TT where {TT<:Real}
    iH1 = findfirst(==(:H1), species)
    return 3.8 * 10^(1) * (1 + xa[iH1]) * exp(lnρ) * (exp(lnT))^(-7/2) + 0.2 * (1 + xa[iH1])  # in cm^2/g
    # return 0.2 * (1 + xa[iH1])
end