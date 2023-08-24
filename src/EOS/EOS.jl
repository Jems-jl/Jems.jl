module EOS

using ..Chem, ..Constants

export AbstractEOS, EOSResults

abstract type AbstractEOS end

@kwdef mutable struct EOSResults{T1}
    T::T1 = 0 
    P::T1 = 0 
    ρ::T1 = 0 
    Prad::T1 = 0 
    μ::T1 = 0 
    α::T1 = 0 
    β::T1 = 0 
    δ::T1 = 0 
    u::T1 = 0 
    cₚ::T1 = 0 
    ∇ₐ::T1 = 0 
    Γ₁::T1 = 0 
end


include("IdealEOS.jl")

end # module EOS
