export set_turb_results!

"""
    struct AdiabaticMLT <: AbstractTurb

Definitions used for MLT without inclusion of radiative losses
"""
struct AdiabaticMLT{T1} <: AbstractTurb where {T1}
    α_MLT::T1
end

set_turb_results!(r::TurbResults, turb::AdiabaticMLT, κ::T1, L::T1, ρ::T1, P::T1, T::T1, r::T1,
                  δ::T1, cₚ::T1, ∇ₐ::T1, m::T2) where {T1<:Real, T2<:Real}
    g = CGRAV*m/r^2
    r.Γ = 0
    r.∇ᵣ = 3κ * L * P / (16π * CRAD * CLIGHT * CGRAV * m * T)
    r.Hₚ = P/(ρ*g)
    if r.∇ᵣ < ∇ₐ
        r.∇ = r.∇_r
        r.v_turb = 0
        r.D_turb = 0
        return
    end
    r.∇ = L/(4πr^2)/cbrt(ρ*cₚ*T*sqrt(g*δ)*(turb.α_MLT*r.Hₚ)^2/(4*sqrt(2)))^2/r.Hₚ + ∇ₐ
    r.v_turb = sqrt(g*δ/Hₚ*(r.∇-r.∇ₐ))*(turb.α_MLT*Hₚ)/(2*sqrt(2))
    r.D_turb = 1/3*r.v_turb*(turb.α_MLT*r.Hₚ)
end