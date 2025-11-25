module Turbulence

using ..Constants

export AbstractTurb, TurbResults

"""
    abstract type AbstractTurb

Abstract supertype from which all turbulence models must derive, _ie_:

`struct MyTurb <: AbstractTurb`
"""
abstract type AbstractTurb end

"""
    mutable struct TurbResults{T1}

Structure that holds various results from the evaluation of the turbulence model of a certain cell.
"""
@kwdef mutable struct TurbResults{T1}
    ∇::T1 = 0  # Actual temperature gradient d lnT/ d lnP
    ∇ᵣ::T1 = 0  # Radiative temperature gradient
    v_turb::T1 = 0  # velocity of turbulent motion
    D_turb::T1 = 0  # mixing coefficient from turbulent motion
    Γ::T1 = 0  # convective efficiency
    Hₚ::T1 = 0  # pressure scale height
end

include("KippenhahnMLT.jl")
include("CoxGiuliMLT.jl")
end