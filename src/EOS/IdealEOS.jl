export set_EOS_resultsTP!

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
    get_EOS_resultsTP(eos, lnT, lnP, xa, species)

computes thermodynamical quantities of a mixture `xa` at temperature `lnT` and pressure `lnP`, given the ideal equation
of state `eos`, list of `species`.
"""
function set_EOS_resultsTP!(eos::IdealEOS, r::EOSResults{TT}, lnT::TT, lnP::TT, xa::AbstractVector{TT},
                           species::Vector{Symbol}) where {TT<:Real}
    # See section 13.2 of Kipp
    r.T = exp(lnT)
    r.P = exp(lnP)
    r.μ = get_μ_IdealEOS(xa, species)
    r.Prad = CRAD * r.T^4 / 3
    r.ρ = r.μ / (CGAS * r.T) * (r.P - r.Prad)

    if eos.include_radiation
        r.β = 1 - r.Prad / r.P  # gas pressure fraction
        if r.β < 0
            throw(DomainError(r.β, "Input values produce negative beta"))
        end
    else
        r.β = 1
    end
    r.α = 1 / r.β
    r.δ = (4 - 3 * r.β) / r.β

    r.u = CGAS * r.T / r.μ * (3 / 2 + 3 * (1 - r.β) / r.β)  # internal energy
    r.cₚ = CGAS / r.μ * (3 / 2 + 3 * (4 + r.β) * (1 - r.β) / r.β^2 + r.δ * r.α)  # heat capacity at constant P
    r.∇ₐ = CGAS * r.δ / (r.β * r.μ * r.cₚ)  # adiabatic temperature gradient
    # first adiabatic exponent dlnP/dlnρ
    if eos.include_radiation
        r.Γ₁ = 1 / (r.α - r.δ * r.∇ₐ)
    else
        r.Γ₁ = 1 / (1 - r.∇ₐ)
    end
    return
end
