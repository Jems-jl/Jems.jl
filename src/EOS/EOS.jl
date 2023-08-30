module EOS

using ..Chem, ..Constants

export AbstractEOS, EOSResults

abstract type AbstractEOS end

@kwdef mutable struct EOSResults{T1<:Real}
    T::T1 = 0  # Temperature (K)
    P::T1 = 0  # Pressure (dyn)
    ρ::T1 = 0  # density (g cm^-3)
    Prad::T1 = 0  # Radiation pressure (dyn)
    μ::T1 = 0  # molecular weight (dimless)
    α::T1 = 0  # == 1/β (dimless)
    β::T1 = 0  # gas pressure fraction; 1 - Prad/P (dimless)
    δ::T1 = 0  # == (4 - 3β) / β (dimless)
    u::T1 = 0  # specific internal energy (erg g^-1)
    cₚ::T1 = 0  # specific heat capacity at constant P (erg g^-1 K^-1)
    ∇ₐ::T1 = 0  # adiabatic temperature gradient dlnT/dlnP (dimless)
    Γ₁::T1 = 0  # first adiabatic exponent dlnP/dlnρ (dimless)
end

include("IdealEOS.jl")

end
