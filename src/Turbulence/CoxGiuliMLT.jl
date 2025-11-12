export set_turb_results!

"""
    struct cgMLT <: AbstractTurb 
Definitions used for MLT with inclusion of radiative losses
"""
struct cgMLT{T1} <: AbstractTurb where {T1}
    α_MLT::T1
    a₀::T1 # geometric factor (9/4 used in Cox & Giuli))
end

function set_turb_results!(turb::cgMLT, res::TurbResults, κ::T1, L::T1, ρ::T1, P::T1, T::T1, r::T1,
                  δ::T1, cₚ::T1, ∇ₐ::T1, m::T2) where {T1<:Real, T2<:Real}
    
    g = CGRAV * m / r^2
    res.Γ = 0 #intitalize the Convective efficiency 
    res.∇ᵣ = 3 * κ * L * P / (16π * CRAD * CLIGHT * CGRAV * m * T^4)
    res.Hₚ = P / (ρ * g)  
    # Convective efficiency factor
    A = sqrt(δ) * cₚ * κ * g * ρ^(5/2) * (turb.α_MLT * res.Hₚ)^2 / (12 * sqrt(2) * CRAD * CLIGHT * T^3 * sqrt(P))
    
    """
    The three parametric equations to be solved for ∇, ∇' and Γ:

        Γ = A * sqrt(∇ - ∇')
        res.∇ᵣ - ∇ = a₀ * A * (∇ - ∇')^(3/2)
        Γ = (∇ - ∇') / (∇' - ∇ₐ)
    """
    
    if res.∇ᵣ < ∇ₐ
        res.∇ = res.∇ᵣ
        res.v_turb = 0
        res.D_turb = 0.0
        return
    end  

    B = (A^2 * (res.∇ᵣ - ∇ₐ) / turb.a₀)^(1/3) 

    # coefficients for the cubic
    a = turb.a₀ * B^2
    b = B
    c = 1
    d = - turb.a₀ * B^2

    Δ₀ = b^2 - 3*a*c
    Δ₁ = 2*b^3 - 9*a*b*c + 27*a^2*d

    # solve cubic
    if Δ₁ < 0
        C = cbrt((Δ₁ - sqrt(Δ₁^2 - 4*Δ₀^3)) / 2)
    else
        C = cbrt((Δ₁ + sqrt(Δ₁^2 - 4*Δ₀^3)) / 2)
    end

    ξ = -1 / (3*a) * (b + C + Δ₀ / C)

    """
    Original equation: ζ^(1/3) + B ζ^(1/2) + a₀ B^2 ζ - a₀ B^2 = 0
    Using ζ = ξ^3 to convert to cubic in ξ
    """
    ζ = ξ^3
    res.∇ = (1 - ζ) * res.∇ᵣ + ζ * ∇ₐ
    res.Γ = B * ξ

    # Turbulent velocity and diffusion
    res.v_turb = 1/(2*sqrt(2)) * turb.α_MLT * (P*δ/ρ)^(1/2) * res.Γ* (1/A)
    res.D_turb = (1/3) * res.v_turb * (turb.α_MLT * res.Hₚ)

end
