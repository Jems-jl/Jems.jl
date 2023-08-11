export get_EOS_resultsTP

"""
    struct IdealEOS <: AbstractEOS

Interface of an Ideal gas equation of state
"""
struct IdealEOS <: AbstractEOS
    include_radiation::Bool
end

"""
    get_μ_IdealEOS(isotope_data, xa, species)

computes the molecular weight of the mixture `xa`, given the `isotope_data` and list of `species`.
"""
function get_μ_IdealEOS(isotope_data::Dict{Symbol,Isotope}, xa::Vector{<:TT},
                        species::Vector{Symbol})::TT where {TT<:Real}
    # See section 4.2 of Kipp
    μ::TT = 0
    for i in eachindex(species)
        # assumes complete ionization!
        μ += xa[i] * (1 + isotope_data[species[i]].Z) / isotope_data[species[i]].mass
    end
    return 1 / μ
end

"""
    get_EOS_resultsTP(eos, isotope_data, lnT, lnP, xa, species)

computes thermodynamical quantities of a mixture `xa` at temperature `lnT` and pressure `lnP`, given the ideal equation
of state `eos`, `isotope_data` and list of `species`.
"""
function get_EOS_resultsTP(eos::IdealEOS, isotope_data::Dict{Symbol,Isotope}, lnT::TT, lnP::TT, xa::Vector{<:TT},
                           species::Vector{Symbol})::Vector{<:TT} where {TT<:Real}
    # See section 13.2 of Kipp
    β::TT = 1
    T = exp(lnT)
    P = exp(lnP)
    μ = get_μ_IdealEOS(isotope_data, xa, species)
    Prad = CRAD * T^4 / 3
    ρ = μ / (CGAS * T) * (P - Prad)

    if eos.include_radiation
        β = 1 - Prad / P  # gas pressure fraction
        if β < 0
            throw(DomainError(β, "Input values produce negative beta"))
        end
    end
    α = 1 / β
    δ = (4 - 3 * β) / β

    u = CGAS * T / μ * (3 / 2 + 3 * (1 - β) / β)  # internal energy
    cₚ = CGAS / μ * (3 / 2 + 3 * (4 + β) * (1 - β) / β^2 + δ * α)  # heat capacity at constant P
    ∇ₐ = CGAS * δ / (β * μ * cₚ)  # adiabatic temperature gradient
    Γ₁::TT = 1  # first adiabatic exponent dlnP/dlnρ
    if eos.include_radiation
        Γ₁ = 1 / (α - δ * ∇ₐ)
    else
        Γ₁ = 1 / (1 - ∇ₐ)
    end

    return [ρ, μ, β, u, cₚ, δ, ∇ₐ, Γ₁]
end
