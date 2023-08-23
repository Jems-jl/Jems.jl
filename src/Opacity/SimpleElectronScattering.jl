export get_opacity_resultsTP

struct SimpleElectronScatteringOpacity <: AbstractOpacity end

"""
    get_opacity_resultsTP(opacity::SimpleElectronScatteringOpacity, lnT::TT, lnP::TT, xa::Vector{<:TT},
                            species::Vector{Symbol})::TT where {TT<:Real}

Evaluates the opacity of the current mixture with mass fractions `xa`, species symbols `species` (both these should be
of length `nspecies`), the natural log of temperature and pressure `lnT`, `lnP`, and the opacity law `opacity`.
"""
function get_opacity_resultsTP(opacity::SimpleElectronScatteringOpacity, lnT::TT, lnP::TT, xa::Vector{<:TT},
                               species::Vector{Symbol})::TT where {TT<:Real}
    iH1 = findfirst(==(:H1), species)
    return 0.2 * (1 + xa[iH1])  # in cm^2/g
end
