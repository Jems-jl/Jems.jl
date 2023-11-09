export set_EOS_resultsTρ!

"""
    struct IdealEOS <: AbstractEOS

Interface of an Ideal gas equation of state
"""
struct IdealEOS <: AbstractEOS
    include_radiation::Bool
end

"""
    get_μ_IdealEOS(xa, species)

computes the molecular weight of the mixture `xa`, given and list of `species`.
"""
function get_μ_IdealEOS(xa::AbstractVector{TT}, species::Vector{Symbol})::TT where {TT<:Real}
    # See section 4.2 of Kipp
    μ::TT = 0
    for i in eachindex(species)
        # assumes complete ionization!
        μ += xa[i] * (1 + Chem.isotope_list[species[i]].Z) / Chem.isotope_list[species[i]].mass
    end
    return 1 / μ
end

"""

    set_EOS_resultsTP!(eos::IdealEOS, r::EOSResults{TT}, lnT::TT, lnP::TT, xa::AbstractVector{TT},
                           species::Vector{Symbol}) where {TT<:Real}

Computes thermodynamical quantities of a mixture `xa` at temperature `lnT` and pressure `lnP`, given the ideal equation
of state `eos` and list of `species`. The results are stored in the EOSResults object `r`.
"""
function set_EOS_resultsTρ!(eos::IdealEOS, r::EOSResults{TT}, lnT::TT, lnρ::TT,
                           xa::AbstractVector{TT}, species::Vector{Symbol}) where {TT<:Real}
    # See section 13.2 of Kipp
    r.lnT = lnT
    r.lnρ = lnρ
    r.T = exp(lnT)
    r.ρ = exp(lnρ)
    r.μ = get_μ_IdealEOS(xa, species)
    if eos.include_radiation
        r.Prad = CRAD * r.T^4 / 3
    else
        r.Prad = 0
    end
    r.P = CGAS * r.T * r.ρ / r.μ + r.Prad
    r.lnP = log(r.P)
    r.β = 1 - r.Prad / r.P  # gas pressure fraction
    r.α = 1 / r.β
    r.δ = (4 - 3 * r.β) / r.β
    r.χ_ρ = r.β
    r.χ_T = 4 - 3 * r.β
    r.u = CGAS * r.T / r.μ * (3 / 2 + 3 * (1 - r.β) / r.β)  # specific internal energy
    r.cₚ = CGAS / r.μ * (3 / 2 + 3 * (4 + r.β) * (1 - r.β) * r.α^2 + r.δ * r.α)  # specific heat capacity at constant P
    r.∇ₐ = CGAS * r.δ / (r.β * r.μ * r.cₚ)  # adiabatic temperature gradient
    r.Γ₁ = 1 / (r.α - r.δ * r.∇ₐ)  # first adiabatic exponent dlnP/dlnρ
end
