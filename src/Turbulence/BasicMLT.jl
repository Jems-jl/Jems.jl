export set_turb_results!

"""
    struct AdiabaticMLT <: AbstractTurb

Definitions used for MLT without inclusion of radiative losses
"""
struct BasicMLT{T1} <: AbstractTurb where {T1}
    α_MLT::T1
end

function set_turb_results!(turb::BasicMLT, res::TurbResults, κ::T1, L::T1, ρ::T1, P::T1, T::T1, r::T1,
                  δ::T1, cₚ::T1, ∇ₐ::T1, m::T2) where {T1<:Real, T2<:Real}
    g = CGRAV*m/r^2
    res.Γ = 0
    res.∇ᵣ = 3κ * L * P / (16π * CRAD * CLIGHT * CGRAV * m * T^4)
    res.Hₚ = P/(ρ*g)
    if res.∇ᵣ < ∇ₐ
        res.∇ = res.∇ᵣ
        res.v_turb = 0
        res.D_turb = 0e0
        return
    end

    U = 3*CRAD*CLIGHT*T^3/(cₚ*ρ^2*κ*(turb.α_MLT*res.Hₚ)^2)*sqrt(8*res.Hₚ/(g*δ))
    W = res.∇ᵣ-∇ₐ

    # coefficients for the cubic (see Kippenhahn equation 7.18)
    a = 1
    b = (8/9 - 3)*U
    c = 3*U^2
    d = -8*U/9*(U^2+W)-U^3

    Δ₀ = b^2 - 3*a*c
    Δ₁ = 2*b^3 -9*a*b*c + 27*a^2*d
    # We choose the sign of the sqrt such that we ensure the number inside the cbrt
    # is not zero. This can mostly happen due to numerical rounding errors.
    if (Δ₁ < 0)
        C = cbrt((Δ₁-sqrt(Δ₁^2-4*Δ₀^3))/2)
    else
        C = cbrt((Δ₁+sqrt(Δ₁^2-4*Δ₀^3))/2)
    end

    ξ = -1/(3*a)*(b+C+Δ₀/C)

    #get ∇ from ξ, equation (7.17) of Kippenhahn
    res.∇ = ∇ₐ + ξ^2 - U^2

    res.v_turb = sqrt(g*δ/(8*res.Hₚ))*(ξ-U)*(turb.α_MLT*res.Hₚ) # 7.6 of Kipp
    res.D_turb = 0* 1/3*res.v_turb*(turb.α_MLT*res.Hₚ)
end