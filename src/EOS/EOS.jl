module EOS

using ..Chem, ..Constants

export AbstractEOS, EOSResults

"""
    abstract type AbstractEOS

Abstract supertype from which all equations of state definitions must derive, _ie_:

`struct MyEOS <: AbstractEOS`
"""
abstract type AbstractEOS end

"""
    mutable struct EOSResults{T1<:Real}

Structure that holds various results from the evaluation of the EOS of a certain cell.
"""
@kwdef mutable struct EOSResults{T1<:Real}
    T::T1 = 0  # Temperature (K)
    P::T1 = 0  # Pressure (dyn)
    ρ::T1 = 0  # density (g cm^-3)
    lnT::T1 = 0
    lnP::T1 = 0
    lnρ::T1 = 0
    Prad::T1 = 0  # Radiation pressure (dyn)
    μ::T1 = 0  # molecular weight (dimless)
    α::T1 = 0  # dlnρ/dlnP at constant temperature (dimless)
    β::T1 = 0  # gas pressure fraction; 1 - Prad/P (dimless)
    δ::T1 = 0  # dlnρ/dlnT at constant pressure (dimless)
    χ_ρ::T1 = 0 # dlnP/dlnρ at constant T.
    χ_T::T1 = 0 # dlnP/dlnT at constant ρ.
    u::T1 = 0  # specific internal energy (erg g^-1)
    cₚ::T1 = 0  # specific heat capacity at constant P (erg g^-1 K^-1)
    ∇ₐ::T1 = 0  # adiabatic temperature gradient; dlnT/dlnP at constant entropy (dimless)
    Γ₁::T1 = 0  # first adiabatic exponent dlnP/dlnρ (dimless)
end

include("IdealEOS.jl")

end
