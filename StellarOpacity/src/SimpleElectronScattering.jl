export get_opacity_resultsTP

struct SimpleElectronScatteringOpacity <: AbstractOpacity

end

function get_opacity_resultsTP(opacity::SimpleElectronScatteringOpacity, isotope_data::Dict{Symbol, Isotope},
        lnT::TT, lnP::TT, xa::Vector{<:TT},species::Vector{Symbol})::TT where{TT<:Real}
    iH1 = findfirst(==(:H1),species)
    return 0.2*(1+xa[iH1]) # in cm^2/g
end
