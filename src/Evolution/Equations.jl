using Jems.DualSupport

"""
    equationHSE(sm::StellarModel, k::Int,
                # varm1::Matrix{TT}, var00::Matrix{TT}, varp1::Matrix{TT},
                # eosm1::EOSResults{TT}, eos00::EOSResults{TT}, eosp1::EOSResults{TT},
                # rates::Matrix{TT},
                # κm1::TT, κ00::TT, κp1::TT)::TT where {TT<:Real}

Default equation of hydrostatic equilibrium. Evaluates for cell `k` of StellarModel `sm` to what degree hydrostatic
equilibrium is satisfied.

# Arguments

  - `sm`: Stellar Model
  - `k`: cell number to consider
#   - `varm1`: Matrix holding the dual numbers of the previous cell (`k-1`)
#   - `var00`: Matrix holding the dual numbers of this cell (`k`)
#   - `varp1`: Matrix holding the dual numbers of the next cell (`k+1`)
#   - `eosm1`: EOSResults object holding the results of the EOS evaluation of the previous cell (`k-1`)
#   - `eos00`: EOSResults object holding the results of the EOS evaluation of the current cell (`k`)
#   - `eosp1`: EOSResults object holding the results of the EOS evaluation of the next cell (`k+1`)
#   - `κm1`: Opacity evaluated at the previous cell (`k-1`)
#   - `κ00`: Opacity evaluated at the current cell (`k`)
#   - `κp1`: Opacity evaluated at the next cell (`k+1`)

# Returns

Residual of comparing dlnP/dm with -GM/4πr^4, where the latter is evaluated at the face of cell `k` and `k+1`.
"""
function equationHSE(sm::StellarModel, k::Int)
    if k == sm.props.nz  # atmosphere boundary condition
        lnP₀ = get_00_dual(sm.props.eos_res[k].lnP)
        r₀ = exp(get_00_dual(sm.props.lnr[k]))
        g₀ = CGRAV * sm.props.mstar / r₀^2
        κ00 = get_00_dual(sm.props.κ[k])
        return lnP₀ - log(2g₀ / (3κ00))  # Eddington gray, ignoring radiation pressure term
    end
    if k==1
        P₀ = get_00_dual(sm.props.eos_res[k].P)
        P₊ = get_p1_dual(sm.props.eos_res[k+1].P)
        ρ₀ = get_00_dual(sm.props.eos_res[k].ρ)
        r₀ = exp(get_00_dual(sm.props.lnr[k]))
        r₊ = exp(get_p1_dual(sm.props.lnr[k+1]))
        rmid₊ = 0.5*(r₀ + r₊)
        return 1-(P₊ + CGRAV*ρ₀^2*(2π/3)*rmid₊^2)/P₀
    end

    #log pressure at cell center of cell k
    lnP₀ = get_00_dual(sm.props.eos_res[k].lnP)
    #log pressure at cell center of cell k+1
    lnP₊ = get_p1_dual(sm.props.eos_res[k+1].lnP)

    lnPface = get_00_dual(sm.props.lnP_face[k])

    r₀ = exp(get_00_dual(sm.props.lnr[k]))
    dm = 0.5*(sm.props.dm[k + 1] + sm.props.dm[k])

    return (exp(lnPface) * (lnP₊ - lnP₀) / dm + CGRAV * sm.props.m[k] / (4π * r₀^4)) /
           (CGRAV * sm.props.m[k] / (4π * r₀^4))
end

"""
    equationT(sm::StellarModel, k::Int,
              varm1::Matrix{TT}, var00::Matrix{TT}, varp1::Matrix{TT},
              eosm1::EOSResults{TT}, eos00::EOSResults{TT}, eosp1::EOSResults{TT},
              rates::Matrix{TT},
              κm1::TT, κ00::TT, κp1::TT)::TT where {TT<:Real}

Default equation of energy transport, evaluated for cell `k` of StellarModel `sm`.

# Arguments

Identical to [`equationHSE`](@ref) for compatibility with StellarModels.TypeStableEquation

# Returns

Residual of comparing dlnT/dm with -∇*GMT/4πr^4P, where the latter is evaluated at the face of cell `k` and `k+1`.
"""
function equationT(sm::StellarModel, k::Int)
    lnT₀ = get_00_dual(sm.props.eos_res[k].lnT)
    if k == sm.props.nz  # atmosphere boundary condition
        L₀ = get_00_dual(sm.props.L[k]) * LSUN
        r₀ = exp(get_00_dual(sm.props.lnr[k]))
        return lnT₀ - log(L₀ / (SIGMA_SB * 4π * r₀^2)) / 4  # Eddington gray, ignoring radiation pressure term
    end
    r₀ = exp(get_00_dual(sm.props.lnr[k]))
    lnT₀ = get_00_dual(sm.props.lnT[k])
    lnT₊ = get_p1_dual(sm.props.lnT[k+1])

    Pface = exp(get_00_dual(sm.props.lnP_face[k]))

    #∇ = get_00_dual(sm.props.turb_res[k].∇)

   # Calculating ∇
    ∇ = get_00_dual(sm.props.∇_face[k])
    L = get_00_dual(sm.props.L[k]) * LSUN
    γ₀ = get_00_dual(sm.props.gamma_turb[k])
    ω = exp(γ₀)
    ρ_face = exp(get_00_dual(sm.props.lnρ_face[k]))
    P_face = exp(get_00_dual(sm.props.lnP_face[k]))
    r₀ = exp(get_00_dual(sm.props.lnr[k]))
    T_face = exp(get_00_dual(sm.props.lnT_face[k]))
    ∇ₐ = get_00_dual(sm.props.∇ₐ_face[k])
    # ∇ᵣ = get_00_dual(sm.props.∇ᵣ_face[k])
    cₚ =  get_00_dual(sm.props.cₚ_face[k])
    κ = get_00_dual(sm.props.κ_face[k])
    m₀ = sm.props.m[k]
    Hₚ = P_face / (ρ_face * CGRAV * m₀ / r₀^2) #defined at face 
    Λ = 1/(1/Hₚ + 1/r₀) 
    k_rad = 16 * SIGMA_SB * T_face^3 / (3 * κ * ρ_face)
    α₂ = ρ_face*cₚ*0.5*sqrt(2/3)*Λ*sqrt(ω)
    ∇ᵣ = 3 * κ * L * P_face / (16π * CRAD * CLIGHT * CGRAV * m₀ * T_face^4)
    SA = (∇ᵣ - ∇ₐ)*(1 + α₂/k_rad)^(-1)
    ∇ = ∇ₐ + SA 
    
#     open("SA_values.txt", "w") do io
#     for k in 1:sm.props.nz
#         println(io, "k=$k, SA=$(SA.value), ∇ = $(∇.value)")
#     end
# end

    
    
    dm = 0.5*(sm.props.dm[k + 1] + sm.props.dm[k])
    
    return ((lnT₊ - lnT₀) / dm + CGRAV * sm.props.m[k] / (4π * r₀^4 * Pface) * ∇) /
        (CGRAV * sm.props.m[k] / (4π * r₀^4 * Pface))

end

"""
    equationLuminosity(sm::StellarModel, k::Int,
                       varm1::Matrix{TT}, var00::Matrix{TT}, varp1::Matrix{TT},
                       eosm1::EOSResults{TT}, eos00::EOSResults{TT}, eosp1::EOSResults{TT},
                       rates::Matrix{TT},
                       κm1::TT, κ00::TT, κp1::TT)::TT where {TT<:Real}

Default equation of energy generation, evaluated for cell `k` of StellarModel `sm`.

# Arguments

Identical to [`equationHSE`](@ref) for compatibility with StellarModels.TypeStableEquation

# Returns

Residual of comparing dL/dm with ϵnuc - cₚ * dT/dt - (δ / ρ) * dP/dt
"""
function equationLuminosity(sm::StellarModel, k::Int)
    L₀ = get_00_dual(sm.props.L[k]) * LSUN
    ρ₀ = get_00_dual(sm.props.eos_res[k].ρ)
    cₚ = get_00_dual(sm.props.eos_res[k].cₚ)
    δ = get_00_dual(sm.props.eos_res[k].δ)
    T₀ = get_00_dual(sm.props.eos_res[k].T)
    P₀ = get_00_dual(sm.props.eos_res[k].P)
    dTdt = (T₀ - get_value(sm.start_step_props.eos_res[k].T)) / sm.props.dt
    dPdt = (P₀ - get_value(sm.start_step_props.eos_res[k].P)) / sm.props.dt

    ϵnuc::typeof(L₀) = 0
    for i in eachindex(sm.network.reactions)
        ϵnuc += get_00_dual(sm.props.rates[k,i])*sm.network.reactions[i].Qvalue
    end
    if k > 1
        L₋ = get_m1_dual(sm.props.L[k-1]) * LSUN
        return ((L₀ - L₋) / sm.props.dm[k] - ϵnuc + cₚ * dTdt - (δ / ρ₀) * dPdt)/(L₀/sm.props.m[k])  # no neutrinos
    else
        return (L₀ / sm.props.dm[k] - ϵnuc + cₚ * dTdt - (δ / ρ₀) * dPdt)/(L₀/sm.props.m[k])  # no neutrinos
    end
end
"""
    equationContinuity(sm::StellarModel, k::Int,
                       varm1::Matrix{TT}, var00::Matrix{TT}, varp1::Matrix{TT},
                       eosm1::EOSResults{TT}, eos00::EOSResults{TT}, eosp1::EOSResults{TT},
                       rates::Matrix{TT},
                       κm1::TT, κ00::TT, κp1::TT)::TT where {TT<:Real}

Default equation of mass continuity, evaluated for cell `k` of StellarModel `sm`.

# Arguments

Identical to [`equationHSE`](@ref) for compatibility with StellarModels.TypeStableEquation

# Returns

Residual of comparing dr^3/dm with 3/(4πρ)
"""
function equationContinuity(sm::StellarModel, k::Int)
    ρ₀ = get_00_dual(sm.props.eos_res[k].ρ)
    r₀ = exp(get_00_dual(sm.props.lnr[k]))
    expected_dr³_dm = 3 / (4π * ρ₀)

    if k > 1  # get inner radius
        r₋ = exp(get_m1_dual(sm.props.lnr[k-1]))
        actual_dr³_dm = (r₀^3 - r₋^3) / sm.props.dm[k]
    else
        actual_dr³_dm = (r₀^3) / sm.props.dm[k]
    end

    return (expected_dr³_dm - actual_dr³_dm) * ρ₀  # times ρ to make eq. dim-less
end

"""
    equation_composition(sm::StellarModel, k::Int, iso_name::Symbol,
                         varm1::Matrix{TT}, var00::Matrix{TT}, varp1::Matrix{TT},
                         eosm1::EOSResults{TT}, eos00::EOSResults{TT}, eosp1::EOSResults{TT},
                         rates::Matrix{TT},
                         κm1::TT, κ00::TT, κp1::TT)::TT where {TT<:Real}

Default equation for composition evolution for isotope `iso_name`, evaluated for cell `k` of StellarModel `sm`.

# Arguments

Identical to [`equationHSE`](@ref) for compatibility with StellarModels.TypeStableEquation

# Returns

Residual of comparing dX_i/dt with its computed reaction rate
"""
function equation_composition(sm::StellarModel, k::Int, iso_name::Symbol)
    # Get mass fraction for this iso
    X00 = get_00_dual(sm.props.xa[k, sm.network.xa_index[iso_name]])
    # We do not use a dual for Xi as this quantity is fixed through the timestep.
    Xi00 = get_value(sm.start_step_props.xa[k, sm.network.xa_index[iso_name]])

    dXdt_nuc::typeof(X00) = 0
    reactions_in = sm.network.species_reactions_in[sm.network.xa_index[iso_name]]
    for reaction_in in reactions_in
        rate = get_00_dual(sm.props.rates[k,reaction_in[1]])
        dXdt_nuc = dXdt_nuc - rate*reaction_in[2]
    end
    reactions_out = sm.network.species_reactions_out[sm.network.xa_index[iso_name]]
    for reaction_out in reactions_out
        rate = get_00_dual(sm.props.rates[k,reaction_out[1]])
        dXdt_nuc = dXdt_nuc + rate*reaction_out[2]
    end
    dXdt_nuc = dXdt_nuc*Chem.isotope_list[iso_name].A*AMU

    #mixing terms
    flux_down::typeof(X00) = 0
    flux_up::typeof(X00) = 0

    # flux term is (4πr^2ρ)^2*D/dm_face, with all quantities evaluated at the face
    # maybe good to experiment with a soft flux limiter that is continuous rather than
    # taking a min
    flux_limiter = sm.opt.physics.flux_limiter
    if k != sm.props.nz 
        Xp1 = get_p1_dual(sm.props.xa[k+1, sm.network.xa_index[iso_name]])
        flux_term_up = get_00_dual(sm.props.flux_term[k])
        flux_term_up = min(flux_limiter*sm.props.dm[k]/sm.props.dt, flux_term_up)
        flux_up = flux_term_up*(Xp1-X00)
    end
    if k != 1
        Xm1 = get_m1_dual(sm.props.xa[k-1, sm.network.xa_index[iso_name]])
        flux_term_down = get_m1_dual(sm.props.flux_term[k-1])
        flux_term_down = min(flux_limiter*sm.props.dm[k]/sm.props.dt, flux_term_down)
        flux_down = flux_term_down*(X00-Xm1)
    end
    dXdt_mix =  (flux_up - flux_down)/(sm.props.dm[k])

    return ((X00 - Xi00) -  (dXdt_nuc + dXdt_mix)*sm.props.dt)
end

# This is unused for now, can be used to implement fractional mixing to a cell that
# is partially convective. Requires :convfrac to be a variable, which represents
# the fraction of a cell that is convective.
function equationConvFrac(sm::StellarModel, k::Int)
    convfrac = get_00_dual(sm.props.convfrac[k])
    if k==1
        diffgrads_00 = get_00_dual(sm.props.turb_res[k].∇ᵣ) - get_00_dual(sm.props.∇ₐ_face[k])
        if diffgrads_00 > 0
            return convfrac - 1.0 # cell is fully convective, so solution is convfrac=1
        else
            return convfrac # cell is fully radiative, so solution is convfrac=0
        end
    end
    if k==sm.props.nz
        diffgrads_m1 = get_m1_dual(sm.props.turb_res[k-1].∇ᵣ) - get_m1_dual(sm.props.∇ₐ_face[k-1])
        if diffgrads_m1 > 0
            return convfrac - 1.0 # cell is fully convective, so solution is convfrac=1
        else
            return convfrac # cell is fully radiative, so solution is convfrac=0
        end
    end
    diffgrads_m1 = get_m1_dual(sm.props.turb_res[k-1].∇ᵣ) - get_m1_dual(sm.props.∇ₐ_face[k-1])
    diffgrads_00 = get_00_dual(sm.props.turb_res[k].∇ᵣ) - get_00_dual(sm.props.∇ₐ_face[k])
    if diffgrads_m1>0 && diffgrads_00>0
        return convfrac - 1.0 # cell is fully convective, so solution is convfrac=1
    elseif diffgrads_m1<0 && diffgrads_00<0
        return convfrac # cell is fully radiative, so solution is convfrac=0
    elseif diffgrads_m1>0
        #convection from the bottom of the cell (diffgrads_00 is negative)
        real_convfrac = diffgrads_m1/(diffgrads_m1-diffgrads_00)
        return convfrac - real_convfrac
    else
        #convection from the top of the cell (diffgrads_m1 is negative)
        real_convfrac = diffgrads_00/(diffgrads_00-diffgrads_m1)
        return convfrac - real_convfrac
    end
end

"""
Equation for turbulent kinetic energy evolution, evaluated for cell `k` of StellarModel `sm`.

Previous function: 
    function equationGammaTurb(sm::StellarModel, k::Int)

    γ = get_00_dual(sm.props.gamma_turb[k])
    ω = exp(γ)
    ∇ₐ = get_00_dual(sm.props.eos_res[k].∇ₐ)
    ∇ = get_00_dual(sm.props.turb_res[k].∇)
    ∇ᵣ = get_00_dual(sm.props.turb_res[k].∇ᵣ)
    cₚ = get_00_dual(sm.props.eos_res[k].cₚ)
    κ = get_00_dual(sm.props.κ[k])
    ρ₀ = get_00_dual(sm.props.eos_res[k].ρ)
    P₀ = get_00_dual(sm.props.eos_res[k].P)
    r₀ = exp(get_00_dual(sm.props.lnr[k]))
    m₀ = sm.props.m[k]
    Hₚ = P₀ / (ρ₀ * CGRAV * m₀ / r₀^2)
    T₀ = get_00_dual(sm.props.eos_res[k].T)
    Λ = 1/(1/Hₚ + 1/r₀)
    τᵣ = cₚ * κ * ρ₀^2 * Λ^2 / (4 * K_BOLTZ * T₀^3)

    α₁ = ∇ₐ * T₀ * Λ * 0.5*sqrt(2/3) * cₚ/ Hₚ^2
    C_d = 8/3 * sqrt(2/3)
    α₂ = ρ₀*cₚ*0.5*sqrt(2/3)*Λ*sqrt(ω)
    SA = (∇ᵣ - ∇ₐ)*(1 + α₂)^(-1)
  
    return  α₁* SA * sqrt(ω) - C_d *ω^(3/2)/Hₚ - ω/τᵣ
end

Since ω = exp(γ) -> dω/dt = ω * dγ/dt, we can rewrite the equation in terms of dγ/dt

"""

# function gammaTurb(sm::StellarModel, k::Int)

#     γ₀ = get_00_dual(sm.props.gamma_turb[k])
#     # ω₀ = exp(γ₀)  
#     ∇ₐ = get_00_dual(sm.props.eos_res[k].∇ₐ)
#     ∇ = get_00_dual(sm.props.turb_res[k].∇)
#     ∇ᵣ = get_00_dual(sm.props.turb_res[k].∇ᵣ)
#     cₚ = get_00_dual(sm.props.eos_res[k].cₚ)
#     κ = get_00_dual(sm.props.κ[k])
#     ρ₀ = get_00_dual(sm.props.eos_res[k].ρ)
#     P₀ = get_00_dual(sm.props.eos_res[k].P)
#     r₀ = exp(get_00_dual(sm.props.lnr[k]))
#     m₀ = sm.props.m[k]
#     Hₚ = P₀ / (ρ₀ * CGRAV * m₀ / r₀^2)
#     T₀ = get_00_dual(sm.props.eos_res[k].T)
#     Λ = 1/(1/Hₚ + 1/r₀)
#     τᵣ = cₚ * κ * ρ₀^2 * Λ^2 / (4 * K_BOLTZ * T₀^3)

#     α₁ = ∇ₐ * T₀ * Λ * 0.5*sqrt(2/3) * cₚ/ Hₚ^2
#     C_d = 8/3 * sqrt(2/3)
#     α₂ = ρ₀*cₚ*0.5*sqrt(2/3)*Λ*sqrt(exp(γ₀))
#     SA = (∇ᵣ - ∇ₐ)*(1 + α₂)^(-1)


    
#     dgammadt = (γ₀ - get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt 

#     # if (k==2)
#     #     @show SA ,sqrt(ω), C_d *ω^(3/2)/Hₚ, ω/τᵣ
#     # end
  
#     return  (dgammadt - α₁* SA/ sqrt(exp(γ₀)) + C_d *sqrt(exp(γ₀))/Λ + 1/τᵣ)/ (exp(γ₀)/sm.props.dt) 
# end

#  function gammaTurb(sm::StellarModel, k::Int)
      
#     γ₀ = get_00_dual(sm.props.gamma_turb[k]) # γ₀ = log(ω) defined at face 
#     ω = exp(γ₀)  # ω defined at face 
#     ω_c = exp(0.5*(get_00_dual(sm.props.gamma_turb[k]) + get_00_dual(sm.props.gamma_turb[k+1]))) # ω defined at cell center
#     r₀ = exp(get_00_dual(sm.props.lnr[k])) # radius of cell face k
#     m₀ = sm.props.m[k] # mass enclosed at cell face k
#     r_c = 0.5 * (exp(get_00_dual(sm.props.lnr[k]))+ exp(get_00_dual(sm.props.lnr[k+1]))) # radius at cell center between k and k+1
#     m_c = 0.5 * (sm.props.m[k] + sm.props.m[k+1]) # mass enclosed at cell center between k and k+1
#     ∇ₐ = get_00_dual(sm.props.∇ₐ_face[k]) #defined at face 
#     ∇ₐ_c = get_00_dual(sm.props.eos_res[k].∇ₐ) # defined at cell center
#     cₚ =  get_00_dual(sm.props.cₚ_face[k]) #defined at face
#     cₚ_c = get_00_dual(sm.props.eos_res[k].cₚ) # defined at cell center
#     κ = get_00_dual(sm.props.κ_face[k]) #defined at face
#     κ_c = get_00_dual(sm.props.κ[k])   # defined at cell center
#     L = get_00_dual(sm.props.L[k]) * LSUN # defined at face
#     L_c = 0.5 * (get_00_dual(sm.props.L[k]) + get_00_dual(sm.props.L[k+1])) * LSUN # defined at cell center

#     ρ₀ = get_00_dual(sm.props.eos_res[k].ρ) # defined at cell center
#     ρ_face = exp(get_00_dual(sm.props.lnρ_face[k])) # defined at face
#     P₀ = get_00_dual(sm.props.eos_res[k].P) # defined at cell center
#     P_face = exp(get_00_dual(sm.props.lnP_face[k])) # defined at face
#     T₀ = get_00_dual(sm.props.eos_res[k].T) # defined at cell center
#     T_face = exp(get_00_dual(sm.props.lnT_face[k])) # defined at face

#     Hₚ = P_face / (ρ_face * CGRAV * m₀ / r₀^2) #defined at face 
#     Hₚ_C = P₀ / (ρ₀ * CGRAV * m_c / r_c^2) # defined at cell center
#     Λ = 1/(1/Hₚ + 1/r₀) #defined at face
#     Λ_c = 1/(1/Hₚ_c + 1/r_c)  # defined at cell center
#     τᵣ = cₚ * κ * ρ_face^2 * Λ^2 / (48 * SIGMA_SB * T_face^3 ) # defined at face
#     τᵣ_c = cₚ_c * κ_c * ρ₀^2 * Λ_c^2 / (48 * SIGMA_SB * T₀^3 ) # defined at cell center
#     c_s  = sqrt(P_face/ρ_face) # defined at face
#     c_s_c  = sqrt(P₀/ρ₀) # defined at cell center
#     k_rad = 16 * SIGMA_SB * T_face^3 / (3 * κ * ρ_face) # defined at face
#     k_rad_c = 16 * SIGMA_SB * T₀^3 / (3 * κ_c * ρ₀) # defined at cell center
#     α₁ = ∇ₐ * T_face * Λ * 0.5*sqrt(2/3) * cₚ/ Hₚ^2 # defined at face
#     α₁_c = ∇ₐ_c * T₀ * Λ_c * 0.5*sqrt(2/3) * cₚ_c/ Hₚ_c^2 # defined at cell center
#     C_d = 8/3 * sqrt(2/3)
#     α₂ = ρ_face*cₚ*0.5*sqrt(2/3)*Λ*sqrt(ω) # defined at face
#     α₂_c = ρ₀*cₚ_c*0.5*sqrt(2/3)*Λ_c*sqrt(ω) # defined at cell center
#     ∇ᵣ = 3 * κ * L * P_face / (16π * CRAD * CLIGHT * CGRAV * m₀ * T_face^4) # defined at face
#     ∇ᵣ_c = 3 * κ_c * L_c * P₀ / (16π * CRAD * CLIGHT * CGRAV * m_c * T₀^4) # defined at cell center
#     SA = (∇ᵣ - ∇ₐ)*(1 + α₂/k_rad)^(-1) # defined at face
#     SA_c = (∇ᵣ_c - ∇ₐ_c)*(1 + α₂_c/k_rad_c)^(-1) # defined at cell center
#     ∇ = ∇ₐ + SA  # defined at face
#     ∇_c = ∇ₐ_c + SA_c # defined at cell center

#     #Diffusion term calculation (Done at cell centres)
    
#     A_i_0 = 4π * α_w * r_c^2 * ρ₀* Λ_c * sqrt(ω_c)
#     A_i_1 







#      if k == sm.props.nz

        
        
#         dgammadt = (γ₀ - get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt
#         # @show  cₚ.value, ∇ₐ.value, ∇ᵣ.value,τᵣ.value
#         return ω*dgammadt - α₁_c* SA_c * sqrt(ω) + C_d *ω^(3/2)/Λ + ω/τᵣ - C_d*(c_s*1e-7)^3/Λ
#     end

#     # for k != nz 
#     Hₚ = P_face / (ρ_face * CGRAV * m₀ / r₀^2) #defined at face 
#     Λ = 1/(1/Hₚ + 1/r₀) 
#     τᵣ = cₚ * κ * ρ_face^2 * Λ^2 / (48 * SIGMA_SB * T_face^3 )
#     c_s  = sqrt(P_face/ρ_face)
#     k_rad = 16 * SIGMA_SB * T_face^3 / (3 * κ * ρ_face)
#     α₁ = ∇ₐ * T_face * Λ * 0.5*sqrt(2/3) * cₚ/ Hₚ^2
#     C_d = 8/3 * sqrt(2/3)
#     α₂ = ρ_face*cₚ*0.5*sqrt(2/3)*Λ*sqrt(ω)
#     ∇ᵣ = 3 * κ * L * P_face / (16π * CRAD * CLIGHT * CGRAV * m₀ * T_face^4)
#     SA = (∇ᵣ - ∇ₐ)*(1 + α₂/k_rad)^(-1)
#     ∇ = ∇ₐ + SA

#     #mixing term 
#     ρ₀ = get_00_dual(sm.props.eos_res[k].ρ)
#     ρ₁ = get_00_dual(sm.props.eos_res[k+1].ρ)
#     P₀ = get_00_dual(sm.props.eos_res[k].P)
#     P₁ = get_00_dual(sm.props.eos_res[k+1].P)
#     r₀ = exp(get_00_dual(sm.props.lnr[k]))
#     r₁ = exp(get_m1_dual(sm.props.lnr[k+1]))
#     # ∇ₐ_c = get_00_dual(sm.props.eos_res[k].∇ₐ)
#     # ∇ᵣ = get_00_dual(sm.props.turb_res[k].∇ᵣ)
#     # cₚ_c = get_00_dual(sm.props.eos_res[k].cₚ)
#     # κ_c = get_00_dual(sm.props.κ[k])
#     Hₚ_c = P₀ / (ρ₀ * CGRAV * (0.5*(sm.props.dm[k] + sm.props.dm[k+1])) / r₀^2)
#     Hₚ_c_1 = P₁ / (ρ₁ * CGRAV * (0.5*(sm.props.dm[k-1] + sm.props.dm[k])) / r₁^2)
#     T₀ = get_00_dual(sm.props.eos_res[k].T)
#     Λ_c= 1/(1/Hₚ_c + 1/r₀)
#     Λ_c_1 = 1/(1/Hₚ_c_1 + 1/r₁) 
#     ω_c_0 = 0.5* (exp(get_00_dual(sm.props.gamma_turb[k])) + exp(get_00_dual(sm.props.gamma_turb[k+1])))
#     ω_c_p = 0.5* (exp(get_00_dual(sm.props.gamma_turb[k-1])) + exp(get_00_dual(sm.props.gamma_turb[k])))
#     ω_c_n = 0.5* (exp(get_00_dual(sm.props.gamma_turb[k+1])) + exp(get_00_dual(sm.props.gamma_turb[k+2])))
#     A_i = (4π * r₀^2 * ρ₀* Λ_c * sqrt(ω_c_0) * ( ω_c_0 - ω_c_p ))/ 0.5*(sm.props.dm[k] + sm.props.dm[k+1])
#     A_i_1 = (4π * r₁^2 * ρ₀* Λ_c_1 * sqrt(ω_c_n) * ( ω_c_n - ω_c_0 ))/ 0.5*(sm.props.dm[k+1] + sm.props.dm[k+2])

#     D_mix = 4π * r₀^2 * (A_i_1 - A_i) / sm.props.dm[k]

#     #∇ = get_00_dual(sm.props.turb_res[k].∇)
#     dgammadt = (γ₀ - get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt
#     # if (k== 1)
#     #     @show (α₁* SA * sqrt(ω)), (ω*dgammadt), (C_d *ω^(3/2)/Hₚ), (ω/τᵣ),  (C_d*(c_s*1e-10)^3/Hₚ)
#     #     @show (ω).value, γ₀.value, ρ_face.value, P_face.value, T_face.value
#     #  end
# #     open("SA_values.txt", "w") do io
# #     for k in 1:sm.props.nz
# #         println(io, "k=$k, SA=$(SA.value), ∇ = $(∇.value), ∇ᵣ=$(∇ᵣ.value), ∇ₐ=$(∇ₐ.value)")
# #     end
# # end

   
#     return  ω*dgammadt - α₁* SA * sqrt(ω) + C_d *ω^(3/2)/Λ + ω/τᵣ - C_d*(c_s*1e-7)^3/Λ + D_mix
# end


"""
function calculate_gamma_turbulance_variables(P, ρ, T, κ, L, r, m_enclosed, ∇ₐ, cₚ, ω)

    C_d = 8 / 3 * sqrt(2 / 3)    
# Pressure Scale Height (Hₚ) and Mixing Length (Λ)
    Hₚ = P / (ρ * CGRAV * m_enclosed / r^2) 
    Λ = 1 / (1 / Hₚ + 1 / r) 
    
    # Radiative Gradient (∇ᵣ)
    ∇ᵣ = (3 * κ * L * P) / (16π * CRAD * CLIGHT * CGRAV * m_enclosed * T^4)
    
    # Radiative Cooling Time Scale (τᵣ)
    τᵣ = (cₚ * κ * ρ^2 * Λ^2) / (48 * SIGMA_SB * T^3)
    
    # Sound Speed (c_s)
    c_s = sqrt(P / ρ)
    
    # Radiative Diffusivity (k_rad)
    k_rad = (16 * SIGMA_SB * T^3) / (3 * κ * ρ)
    
    # MLT Coefficients
    α₂ = ρ * cₚ * 0.5 * sqrt(2 / 3) * Λ * sqrt(ω)
    α₁ = ∇ₐ * T * Λ * 0.5 * sqrt(2 / 3) * cₚ / Hₚ^2
    
    # Super-Adiabatic Excess (SA)
    SA = (∇ᵣ - ∇ₐ) * (1 + α₂ / k_rad)^(-1)
    
    return (; Hₚ, Λ, ∇ᵣ, τᵣ, c_s, k_rad, α₁, α₂, SA, C_d)

end 

function calculate_flux_F(sm::StellarModel, i::Int)

    if i == sm.props.nz + 1 || i == 1
        return 0.0
    end

    P_c_i = get_00_dual(sm.props.eos_res[i].P)
    ρ_c_i = get_00_dual(sm.props.eos_res[i].ρ)
    T_c_i = get_00_dual(sm.props.eos_res[i].T)
    κ_c_i = get_00_dual(sm.props.κ[i])
    ∇ₐ_c_i = get_00_dual(sm.props.eos_res[i].∇ₐ)
    cₚ_c_i = get_00_dual(sm.props.eos_res[i].cₚ)
    L_c_i = 0.5 * (get_00_dual(sm.props.L[i]) + get_00_dual(sm.props.L[i+1])) * LSUN # Avg L
    r_c_i = 0.5 * (exp(get_00_dual(sm.props.lnr[i])) + exp(get_00_dual(sm.props.lnr[i+1]))) # Cell center radius
    # m_c_i = 0.5 * (sm.props.m[i] + sm.props.m[i+1]) # Cell center enclosed mass
    ω_i = exp(get_00_dual(sm.props.gamma_turb[i]))
    ω_i_minus1 = exp(get_00_dual(sm.props.gamma_turb[i-1]))
    # ω_i_plus1 = exp(get_00_dual(sm.props.gamma_turb[i+1]))
    ω_i_center = 0.5 * (ω_i + ω_i_minus1)
    MLT_Set = calculate_gamma_turbulance_variables(P_c_i, ρ_c_i, T_c_i, κ_c_i, L_c_i, r_c_i, m_c_i, ∇ₐ_c_i, cₚ_c_i, ω_i)
    A_i = 4π * α_w * r_c_i^2 * ρ_c_i * MLT_Set.Λ * sqrt(ω_i_center)

    return A_i * (ω_i - ω_i_minus1)/ (0.5 * (sm.props.dm[i] + sm.props.dm[i-1]))
    
end 

function gammaTurb(sm::StellarModel, k::Int)
    
    #FACE TERMS 
    γ₀ = get_00_dual(sm.props.gamma_turb[k]) # gamma at face k
    ω = exp(γ₀)                              # omega at face k
    m_enclosed = sm.props.m[k]               # mass enclosed at face k
    r_face = exp(get_00_dual(sm.props.lnr[k])) # radius at face k
    P_f = exp(get_00_dual(sm.props.lnP_face[k]))
    ρ_f = exp(get_00_dual(sm.props.lnρ_face[k]))
    T_f = exp(get_00_dual(sm.props.lnT_face[k]))
    κ_f = get_00_dual(sm.props.κ_face[k])
    ∇ₐ_f = get_00_dual(sm.props.∇ₐ_face[k])
    cₚ_f = get_00_dual(sm.props.cₚ_face[k])
    L_f = get_00_dual(sm.props.L[k]) * LSUN

    #CENTRE TERMS
    γ_c = 0.5*(get_00_dual(sm.props.gamma_turb[k]) + get_00_dual(sm.props.gamma_turb[k-1])) # gamma at center between k and k-1
    ω_c = exp(γ_c)
    m_c = 0.5 * (sm.props.m[k] + sm.props.m[k-1])
    r_c = 0.5 * (exp(get_00_dual(sm.props.lnr[k])) + exp(get_00_dual(sm.props.lnr[k-1])))
    P_c = get_00_dual(sm.props.eos_res[k].P)
    ρ_c = get_00_dual(sm.props.eos_res[k].ρ)
    T_c = get_00_dual(sm.props.eos_res[k].T)
    κ_c = get_00_dual(sm.props.κ[k])
    ∇ₐ_c = get_00_dual(sm.props.eos_res[k].∇ₐ)
    cₚ_c = get_00_dual(sm.props.eos_res[k].cₚ)
    L_c = 0.5 * (get_00_dual(sm.props.L[k]) + get_00_dual(sm.props.L[k-1])) * LSUN

    


    Face_k_MLT = calculate_gamma_turbulance_variables(P_f, ρ_f, T_f, κ_f, L_f, r_face, m_enclosed, ∇ₐ_f, cₚ_f, ω)
    Centre_k_MLT = calculate_gamma_turbulance_variables(P_c, ρ_c, T_c, κ_c, L_c, r_c, m_c, ∇ₐ_c, cₚ_c, ω_c)


    if k == 1
        
        mixing_term = calculate_flux_F(sm, 2)/(sm.props.dm[1])

        dgammadt = (γ₀ - get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt

        return ω * dgammadt - Face_k_MLT.α₁ * Face_k_MLT.SA * sqrt(ω) + Face_k_MLT.C_d * ω^(3/2) / Face_k_MLT.Λ + ω / Face_k_MLT.τᵣ - Face_k_MLT.C_d * (Face_k_MLT.c_s * 1e-7)^3 / Face_k_MLT.Λ + mixing_term
    end

    if k == sm.props.nz
        mixing_term = -calculate_flux_F(sm, sm.props.nz)/(sm.props.dm[sm.props.nz])

        dgammadt = (γ_c - 0.5*(get_value(sm.start_step_props.gamma_turb[k]) + get_value(sm.start_step_props.gamma_turb[k-1]))) / sm.props.dt

        return ω_c * dgammadt - Centre_k_MLT.α₁ * Centre_k_MLT.SA * sqrt(ω_c) + Centre_k_MLT.C_d * ω_c^(3/2) / Centre_k_MLT.Λ + ω_c / Centre_k_MLT.τᵣ - Centre_k_MLT.C_d * (Centre_k_MLT.c_s * 1e-7)^3 / Centre_k_MLT.Λ + mixing_term
    end

    mixing_term = (calculate_flux_F(sm, k + 1) - calculate_flux_F(sm, k)) / (sm.props.dm[k])
    dgammadt = (γ₀ - get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt

    return ω * dgammadt - Face_k_MLT.α₁ * Face_k_MLT.SA * sqrt(ω) + Face_k_MLT.C_d * ω^(3/2) / Face_k_MLT.Λ + ω / Face_k_MLT.τᵣ - Face_k_MLT.C_d * (Face_k_MLT.c_s * 1e-7)^3 / Face_k_MLT.Λ + mixing_term
end
"""

# # NOTE: Assume α_w is a globally defined constant.
# # NOTE: The physical terms Hₚ_C, m_c, r_c, etc., are usually defined between cell k and k-1.

# --- Helper function remains the same ---
# function calculate_gamma_turbulance_variables(P, ρ, T, κ, L, r, m_enclosed, ∇ₐ, cₚ, ω)
#     C_d = 8 / 3 * sqrt(2 / 3)    
#     Hₚ = P / (ρ * CGRAV * m_enclosed / r^2) 
#     Λ = 1 / (1 / Hₚ + 1 / r) 
#     ∇ᵣ = (3 * κ * L * P) / (16π * CRAD * CLIGHT * CGRAV * m_enclosed * T^4)
#     τᵣ = (cₚ * κ * ρ^2 * Λ^2) / (48 * SIGMA_SB * T^3)
#     c_s = sqrt(P / ρ)
#     k_rad = (16 * SIGMA_SB * T^3) / (3 * κ * ρ)
#     α₂ = ρ * cₚ * 0.5 * sqrt(2 / 3) * Λ * sqrt(ω)
#     α₁ = ∇ₐ * T * Λ * 0.5 * sqrt(2 / 3) * cₚ / Hₚ^2
#     SA = (∇ᵣ - ∇ₐ) * (1 + α₂ / k_rad)^(-1)
#     return (; Hₚ, Λ, ∇ᵣ, τᵣ, c_s, k_rad, α₁, α₂, SA, C_d)
# end 


# function calculate_flux_F(sm::StellarModel, i::Int)

#     if i == sm.props.nz + 1 || i == 1
#         return 0.0
#     end
#     α_w = 0.25
#     P_c_i = get_00_dual(sm.props.eos_res[i].P)
#     ρ_c_i = get_00_dual(sm.props.eos_res[i].ρ)
#     T_c_i = get_00_dual(sm.props.eos_res[i].T)
#     κ_c_i = get_00_dual(sm.props.κ[i])
#     ∇ₐ_c_i = get_00_dual(sm.props.eos_res[i].∇ₐ)
#     cₚ_c_i = get_00_dual(sm.props.eos_res[i].cₚ)
#     L_c_i = 0.5 * (get_00_dual(sm.props.L[i]) + get_m1_dual(sm.props.L[i-1])) * LSUN # Center L using i and i-1
#     r_c_i = 0.5 * (exp(get_00_dual(sm.props.lnr[i])) + exp(get_m1_dual(sm.props.lnr[i-1]))) # Center radius using i and i-1
#     m_c_i = 0.5 * (sm.props.m[i] + sm.props.m[i-1]) # Center enclosed mass
    
#     # Omega variables
#     ω_i = exp(get_00_dual(sm.props.gamma_turb[i]))
#     ω_i_minus1 = exp(get_m1_dual(sm.props.gamma_turb[i-1]))
#     ω_i_center = 0.5 * (ω_i + ω_i_minus1) # Center omega based on cell average

#     MLT_Set = calculate_gamma_turbulance_variables(P_c_i, ρ_c_i, T_c_i, κ_c_i, L_c_i, r_c_i, m_c_i, ∇ₐ_c_i, cₚ_c_i, ω_i_center)
#     A_i = (4π * r_c_i^2)^2 * ρ_c_i * MLT_Set.Λ * sqrt(ω_i_center)* α_w 

    
#     return A_i * (ω_i - ω_i_minus1) / (0.5 * (sm.props.dm[i] + sm.props.dm[i-1]))
    
# end 

# function calculate_flux_F_p1(sm::StellarModel, i::Int)

#     if i == sm.props.nz + 1 || i == 1
#         return 0.0
#     end
#     α_w = 0.25
#     P_c_i = get_00_dual(sm.props.eos_res[i].P)
#     ρ_c_i = get_00_dual(sm.props.eos_res[i].ρ)
#     T_c_i = get_00_dual(sm.props.eos_res[i].T)
#     κ_c_i = get_00_dual(sm.props.κ[i])
#     ∇ₐ_c_i = get_00_dual(sm.props.eos_res[i].∇ₐ)
#     cₚ_c_i = get_00_dual(sm.props.eos_res[i].cₚ)
#     L_c_i = 0.5 * (get_00_dual(sm.props.L[i]) + get_m1_dual(sm.props.L[i-1])) * LSUN # Center L using i and i-1
#     r_c_i = 0.5 * (exp(get_00_dual(sm.props.lnr[i])) + exp(get_m1_dual(sm.props.lnr[i-1]))) # Center radius using i and i-1
#     m_c_i = 0.5 * (sm.props.m[i] + sm.props.m[i-1]) # Center enclosed mass
    
#     # Omega variables
#     ω_i = exp(get_00_dual(sm.props.gamma_turb[i-1]))
#     ω_i_p1 = exp(get_p1_dual(sm.props.gamma_turb[i]))
#     ω_i_center = 0.5 * (ω_i + ω_i_p1) # Center omega based on cell average

#     MLT_Set = calculate_gamma_turbulance_variables(P_c_i, ρ_c_i, T_c_i, κ_c_i, L_c_i, r_c_i, m_c_i, ∇ₐ_c_i, cₚ_c_i, ω_i_center)
#     A_i = (4π * r_c_i^2)^2 * ρ_c_i * MLT_Set.Λ * sqrt(ω_i_center)* α_w 

    
#     return A_i * (ω_i_p1 - ω_i) / (0.5 * (sm.props.dm[i] + sm.props.dm[i-1]))
    
# end 



# function gammaTurb(sm::StellarModel, k::Int)
    
#     # --- A. FACE TERMS (Used for fluid residual in k=1 and k=interior) ---
#     γ₀ = get_00_dual(sm.props.gamma_turb[k]) # gamma at face k
#     ω = exp(γ₀)                              # omega at face k
#     m_enclosed = sm.props.m[k]               # mass enclosed at face k (m₀)
#     r_face = exp(get_00_dual(sm.props.lnr[k])) # radius at face k (r₀)
#     P_f = exp(get_00_dual(sm.props.lnP_face[k]))
#     ρ_f = exp(get_00_dual(sm.props.lnρ_face[k]))
#     T_f = exp(get_00_dual(sm.props.lnT_face[k]))
#     κ_f = get_00_dual(sm.props.κ_face[k])
#     ∇ₐ_f = get_00_dual(sm.props.∇ₐ_face[k])
#     cₚ_f = get_00_dual(sm.props.cₚ_face[k])
#     L_f = get_00_dual(sm.props.L[k]) * LSUN
    
#     Face_k_MLT = calculate_gamma_turbulance_variables(P_f, ρ_f, T_f, κ_f, L_f, r_face, m_enclosed, ∇ₐ_f, cₚ_f, ω)


#     if k == 1
#         # Mixing term: Flux at k+1 (outer face) divided by cell mass dm[1]. Flux at k=1 is zero.
#         mixing_term = calculate_flux_F_p1(sm, 2) / (sm.props.dm[1])
#         dgammadt = (γ₀ - get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt
        
#         residual = (ω * dgammadt - Face_k_MLT.α₁ * Face_k_MLT.SA * sqrt(ω) + 
#                           Face_k_MLT.C_d * ω^(3/2) / Face_k_MLT.Λ + ω / Face_k_MLT.τᵣ - 
#                           Face_k_MLT.C_d * (Face_k_MLT.c_s * 1e-7)^3 / Face_k_MLT.Λ)
        
#         return residual + mixing_term
#     end



#     if k == sm.props.nz
       
#         γ_c = 0.5*(get_00_dual(sm.props.gamma_turb[k]) + get_m1_dual(sm.props.gamma_turb[k-1])) 
#         ω_c = exp(γ_c)
#         m_c = 0.5 * (sm.props.m[k] + sm.props.m[k-1])
#         r_c = 0.5 * (exp(get_00_dual(sm.props.lnr[k])) + exp(get_m1_dual(sm.props.lnr[k-1])))
#         P_c = get_00_dual(sm.props.eos_res[k].P)
#         ρ_c = get_00_dual(sm.props.eos_res[k].ρ)
#         T_c = get_00_dual(sm.props.eos_res[k].T)
#         κ_c = get_00_dual(sm.props.κ[k])
#         ∇ₐ_c = get_00_dual(sm.props.eos_res[k].∇ₐ)
#         cₚ_c = get_00_dual(sm.props.eos_res[k].cₚ)
#         L_c = 0.5 * (get_00_dual(sm.props.L[k]) + get_m1_dual(sm.props.L[k-1])) * LSUN
        
#         Centre_k_MLT = calculate_gamma_turbulance_variables(P_c, ρ_c, T_c, κ_c, L_c, r_c, m_c, ∇ₐ_c, cₚ_c, ω_c)
        
        
#         mixing_term = -calculate_flux_F(sm, sm.props.nz) / (sm.props.dm[sm.props.nz]) 
        
        
#         dgammadt = (γ_c - 0.5*(get_value(sm.start_step_props.gamma_turb[k]) + get_value(sm.start_step_props.gamma_turb[k-1]))) / sm.props.dt
        
#         residual = (ω_c * dgammadt - Centre_k_MLT.α₁ * Centre_k_MLT.SA * sqrt(ω_c) + 
#                           Centre_k_MLT.C_d * ω_c^(3/2) / Centre_k_MLT.Λ + ω_c / Centre_k_MLT.τᵣ - 
#                           Centre_k_MLT.C_d * (Centre_k_MLT.c_s * 1e-7)^3 / Centre_k_MLT.Λ)
                          
#         return residual + mixing_term
#     end

#     # INTERIOR CELLS (1 < k < nz)

#     mixing_term = (calculate_flux_F_p1(sm, k + 1) - calculate_flux_F(sm, k)) / (sm.props.dm[k])
#     dgammadt = (γ₀ - get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt
    
#     residual = (ω * dgammadt - Face_k_MLT.α₁ * Face_k_MLT.SA * sqrt(ω) + 
#                       Face_k_MLT.C_d * ω^(3/2) / Face_k_MLT.Λ + ω / Face_k_MLT.τᵣ - 
#                       Face_k_MLT.C_d * (Face_k_MLT.c_s * 1e-7)^3 / Face_k_MLT.Λ)
                      
#     return residual + mixing_term
# end

# # NOTE: Constants like LSUN, CGRAV, CRAD, CLIGHT, SIGMA_SB, and α_w are assumed to be defined globally.






function gammaTurb(sm::StellarModel, k::Int)

###Constants###
C_d = 8/3 * sqrt(2/3)
α_w = 0.25 


"""
mixing term = 1/ dm(i,face) [ A(i+1) * (ω(i+1)- ω(i))/ dm(i+1,cell) - A(i) * (ω(i)- ω(i-1))/ dm(i,cell)]
A(i) = (4πr^2(c,i))^2 * ρ(i,c) * Λ(i,c)^2 * √ω(i,c) { X(i,c) = 0.5(X(i,f) + X(i-1,f))}
A(1) = A(nz+1)= 0
"""

## Outer boundary condition (k = 1)
if k == 1
    ### face Values required for calculating all the other terms except mixing term ###
    γ_face_00 = get_00_dual(sm.props.gamma_turb[k])
    ω_face_00 = exp(γ_face_00)  ###
    dm_cell_00 = sm.props.dm[k] 
    m_cell_00 = sm.props.m[k]
    r_face_00 = exp(get_00_dual(sm.props.lnr[k]))
    P_face_00 = exp(get_00_dual(sm.props.lnP_face[k]))
    ρ_face_00 = exp(get_00_dual(sm.props.lnρ_face[k]))
    T_face_00 = exp(get_00_dual(sm.props.lnT_face[k]))
    κ_face_00 = get_00_dual(sm.props.κ_face[k])
    ∇ₐ_face_00 = get_00_dual(sm.props.∇ₐ_face[k])
    cₚ_face_00 = get_00_dual(sm.props.cₚ_face[k])
    L_face_00 = get_00_dual(sm.props.L[k]) * LSUN

    Hₚ_face_00 = P_face_00 / (ρ_face_00 * m_cell_00 * CGRAV / r_face_00^2)
    Λ_face_00 = 1 / (1 / Hₚ_face_00 + 1 / r_face_00)
    ∇ᵣ_face_00 = (3 * κ_face_00 * L_face_00 * P_face_00) / (16π * CRAD * CLIGHT * CGRAV * m_cell_00 * T_face_00^4)
    τᵣ_face_00 = (cₚ_face_00 * κ_face_00 * ρ_face_00^2 * Λ_face_00^2) / (48 * SIGMA_SB * T_face_00^3)
    c_s_face_00 = sqrt(P_face_00 / ρ_face_00)
    k_rad_face_00 = (16 * SIGMA_SB * T_face_00^3) / (3 * κ_face_00 * ρ_face_00)  ##
    α₂_face_00 = ρ_face_00 * cₚ_face_00 * 0.5 * sqrt(2 / 3) * Λ_face_00 * sqrt(ω_face_00)
    α₁_face_00 = ∇ₐ_face_00 * T_face_00 * Λ_face_00 * 0.5 * sqrt(2 / 3) * cₚ_face_00 / Hₚ_face_00^2
    SA_face_00 = (∇ᵣ_face_00 - ∇ₐ_face_00) * (1 + α₂_face_00 / k_rad_face_00)^(-1)
    dgammadt_face_00 = (γ_face_00- get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt

    """
    mixing term = 1/ dm(i,face) [ A(2) * (ω(2)- ω(1))/ dm(2,cell) ]

    """

    ### p1 terms needed as F_1 = 0 and F_2 = 0 is the next cell so everything is wrt to i+1 so the p1 cell, also mixing term calcualted at cell centre ###
    P_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].P)
    ρ_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].ρ)
    T_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].T)
    dm_face_p1 = 0.5*(sm.props.dm[k] + sm.props.dm[k+1])##
    dm_cell_p1 = sm.props.dm[k+1] 
    m_face_p1 = 0.5*(sm.props.m[k] + sm.props.m[k+1])
    L_cc_p1 = 0.5*(get_p1_dual(sm.props.L[k+1]) + get_00_dual(sm.props.L[k])) * LSUN
    r_cc_p1 = 0.5*(exp(get_p1_dual(sm.props.lnr[k+1])) + exp(get_00_dual(sm.props.lnr[k])))
    Hₚ_cc_p1 = P_cc_p1 / (ρ_cc_p1 * m_face_p1 * CGRAV / r_cc_p1^2)
    ω_cc_p1 = exp(0.5*(get_p1_dual(sm.props.gamma_turb[k+1]) + get_00_dual(sm.props.gamma_turb[k])))

    Λ_cc_p1 = 1 / (1 / Hₚ_cc_p1 + 1 / r_cc_p1)
    A_p1 = (4π  *ρ_cc_p1* r_cc_p1^2)^2 * Λ_cc_p1 * α_w * sqrt(ω_cc_p1)

    F_p1 = (A_p1 / dm_cell_p1) * (exp(get_p1_dual(sm.props.gamma_turb[k+1])) - exp(get_00_dual(sm.props.gamma_turb[k]))) 

    # Different terms for residual at k = 1
    mixing_term =  (F_p1/  dm_face_p1)  
    omega_var_term = ω_face_00 * dgammadt_face_00  #ω_face_00 *
    source_term = α₁_face_00 * SA_face_00 *sqrt(ω_face_00)
    turb_dissipation_term = C_d * (ω_face_00)^(3/2) / Λ_face_00
    rad_dissipation_term = ω_face_00 / τᵣ_face_00
    excess_term = C_d * (c_s_face_00 * 1e-7)^3 / Λ_face_00

    return omega_var_term - source_term + turb_dissipation_term + rad_dissipation_term - excess_term - mixing_term
end



## Inner boundary condition (k = sm.props.nz)
if k == sm.props.nz

### At the centre, all values are calculated at the cell centre instead of face ###

    γ_cc_00 = get_00_dual(sm.props.gamma_turb[k]) 
    ω_cc_00 = exp(γ_cc_00) ###
    dm_face_00 = 0.5*(sm.props.dm[k] + sm.props.dm[k-1])
    m_face_00 = 0.5*(sm.props.m[k] + sm.props.m[k-1])
    r_cc_00 = exp(get_00_dual(sm.props.lnr[k]))
    P_cc_00 = get_00_dual(sm.props.eos_res[k].P)
    ρ_cc_00 = get_00_dual(sm.props.eos_res[k].ρ)
    T_cc_00 = get_00_dual(sm.props.eos_res[k].T)
    κ_cc_00 = get_00_dual(sm.props.κ[k])
    ∇ₐ_cc_00 = get_00_dual(sm.props.eos_res[k].∇ₐ)
    cₚ_cc_00 = get_00_dual(sm.props.eos_res[k].cₚ)
    L_cc_00 = get_00_dual(sm.props.L[k])* LSUN

    Hₚ_cc_00 = P_cc_00 / (ρ_cc_00 * m_face_00 * CGRAV / r_cc_00^2)
    Λ_cc_00 = 1 / (1 / Hₚ_cc_00 + 1 / r_cc_00)
    ∇ᵣ_cc_00 = (3 * κ_cc_00 * L_cc_00 * P_cc_00) / (16π * CRAD * CLIGHT * CGRAV * m_face_00 * T_cc_00^4)
    τᵣ_cc_00 = (cₚ_cc_00 * κ_cc_00 * ρ_cc_00^2 * Λ_cc_00^2) / (48 * SIGMA_SB * T_cc_00^3)
    c_s_cc_00 = sqrt(P_cc_00 / ρ_cc_00)
    k_rad_cc_00 = (16 * SIGMA_SB * T_cc_00^3) / (3 * κ_cc_00 * ρ_cc_00)  ##
    α₂_cc_00 = ρ_cc_00 * cₚ_cc_00 * 0.5 * sqrt(2 / 3) * Λ_cc_00 * sqrt(ω_cc_00)
    α₁_cc_00 = ∇ₐ_cc_00 * T_cc_00 * Λ_cc_00 * 0.5 * sqrt(2 / 3) * cₚ_cc_00 / Hₚ_cc_00^2
    SA_cc_00 = (∇ᵣ_cc_00 - ∇ₐ_cc_00) * (1 + α₂_cc_00 / k_rad_cc_00)^(-1)    
    dgammadt_cc_00 = (γ_cc_00- (get_value(sm.start_step_props.gamma_turb[k]) )) / sm.props.dt

    ### Since the F_nz+1 = 0, and only F_nz is present which uses the 00 cell ###
    A_00 = (4π * ρ_cc_00 * r_cc_00^2)^2 * Λ_cc_00 * α_w * sqrt(ω_cc_00)
    F_00 = (A_00 / sm.props.dm[k]) * (exp(get_00_dual(sm.props.gamma_turb[k])) - exp(get_m1_dual(sm.props.gamma_turb[k-1])))

    ## Different terms for residual at k = nz 
    mixing_term =  -(F_00/ sm.props.dm[k])
    omega_var_term = ω_cc_00* dgammadt_cc_00 #ω_cc_00*
    source_term = α₁_cc_00 * SA_cc_00 * sqrt(ω_cc_00)
    turb_dissipation_term = C_d * (ω_cc_00)^(3/2) / Λ_cc_00
    rad_dissipation_term = ω_cc_00 / τᵣ_cc_00
    excess_term = C_d * (c_s_cc_00 * 1e-7)^3 / Λ_cc_00   
    # @show k
    # @show (γ_cc_00).value, (ω_cc_00).value
    # @show (dm_face_00), (m_face_00), (r_cc_00).value
    # @show (Hₚ_cc_00).value, (Λ_cc_00).value
    # @show (dgammadt_cc_00).value
    # @show (mixing_term).value, (omega_var_term).value, (source_term).value,
    #     (turb_dissipation_term).value, (rad_dissipation_term).value, (excess_term).value
    return omega_var_term - source_term + turb_dissipation_term + rad_dissipation_term - excess_term - mixing_term


end 



### Other Calculations : 1 < k < nz ###
begin 

    ### face Values are required for all other terms except the mixing term ###
    γ_face_00 = get_00_dual(sm.props.gamma_turb[k])
    ω_face_00 = exp(γ_face_00) ###
    dm_cell_00 = sm.props.dm[k] 
    m_cell_00 = sm.props.m[k]
    r_face_00 = exp(get_00_dual(sm.props.lnr[k]))
    P_face_00 = exp(get_00_dual(sm.props.lnP_face[k]))
    ρ_face_00 = exp(get_00_dual(sm.props.lnρ_face[k]))
    T_face_00 = exp(get_00_dual(sm.props.lnT_face[k]))
    κ_face_00 = get_00_dual(sm.props.κ_face[k])
    ∇ₐ_face_00 = get_00_dual(sm.props.∇ₐ_face[k])
    cₚ_face_00 = get_00_dual(sm.props.cₚ_face[k])
    L_face_00 = get_00_dual(sm.props.L[k]) * LSUN

    Hₚ_face_00 = P_face_00 / (ρ_face_00 * m_cell_00 * CGRAV / r_face_00^2)
    Λ_face_00 = 1 / (1 / Hₚ_face_00 + 1 / r_face_00)
    ∇ᵣ_face_00 = (3 * κ_face_00 * L_face_00 * P_face_00) / (16π * CRAD * CLIGHT * CGRAV * m_cell_00 * T_face_00^4)
    τᵣ_face_00 = (cₚ_face_00 * κ_face_00 * ρ_face_00^2 * Λ_face_00^2) / (48 * SIGMA_SB * T_face_00^3)
    c_s_face_00 = sqrt(P_face_00 / ρ_face_00)
    k_rad_face_00 = (16 * SIGMA_SB * T_face_00^3) / (3 * κ_face_00 * ρ_face_00)  ##
    α₂_face_00 = ρ_face_00 * cₚ_face_00 * 0.5 * sqrt(2 / 3) * Λ_face_00 * sqrt(ω_face_00)
    α₁_face_00 = ∇ₐ_face_00 * T_face_00 * Λ_face_00 * 0.5 * sqrt(2 / 3) * cₚ_face_00 / Hₚ_face_00^2
    SA_face_00 = (∇ᵣ_face_00 - ∇ₐ_face_00) * (1 + α₂_face_00 / k_rad_face_00)^(-1)
    dgammadt_face_00 = (γ_face_00- get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt
    


    ## cell centre values required for mixing term ##
    ## Both 00 and p1 terms are needed for F_i and F_i+1

    ## 00 values ##
    γ_cc_00 = 0.5*(get_00_dual(sm.props.gamma_turb[k]) + get_m1_dual(sm.props.gamma_turb[k-1])) 
    ω_cc_00 = exp(γ_cc_00) ###
    dm_face_00 = 0.5*(sm.props.dm[k] + sm.props.dm[k-1]) ##
    m_face_00 = 0.5*(sm.props.m[k] + sm.props.m[k-1])
    r_cc_00 = 0.5*(exp(get_00_dual(sm.props.lnr[k])) + exp(get_m1_dual(sm.props.lnr[k-1])))
    P_cc_00 = get_00_dual(sm.props.eos_res[k].P)
    ρ_cc_00 = get_00_dual(sm.props.eos_res[k].ρ)
    L_cc_00 = 0.5*(get_00_dual(sm.props.L[k]) + get_m1_dual(sm.props.L[k-1])) * LSUN
    Hₚ_cc_00 = P_cc_00 / (ρ_cc_00 * m_face_00 * CGRAV / r_cc_00^2)
    Λ_cc_00 = 1 / (1 / Hₚ_cc_00 + 1 / r_cc_00)
    A_00 = (4π  *ρ_cc_00 * r_cc_00^2)^2 * Λ_cc_00 * α_w * sqrt(ω_cc_00)
    F_00 = (A_00 / sm.props.dm[k]) * (exp(get_00_dual(sm.props.gamma_turb[k])) - exp(get_m1_dual(sm.props.gamma_turb[k-1])))

    ## P1 values ##
    P_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].P)
    ρ_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].ρ)
    T_cc_p1 = get_p1_dual(sm.props.eos_res[k+1].T)
    dm_face_p1 = 0.5*(sm.props.dm[k] + sm.props.dm[k+1])
    m_face_p1 = 0.5*(sm.props.m[k] + sm.props.m[k+1]) ##
    L_cc_p1 = 0.5*(get_p1_dual(sm.props.L[k+1]) + get_00_dual(sm.props.L[k])) * LSUN
    r_cc_p1 = 0.5*(exp(get_p1_dual(sm.props.lnr[k+1])) + exp(get_00_dual(sm.props.lnr[k])))
    Hₚ_cc_p1 = P_cc_p1 / (ρ_cc_p1 * m_face_p1 * CGRAV / r_cc_p1^2)
    Λ_cc_p1 = 1 / (1 / Hₚ_cc_p1 + 1 / r_cc_p1)
    ω_cc_p1 = exp(0.5*(get_p1_dual(sm.props.gamma_turb[k+1]) + get_00_dual(sm.props.gamma_turb[k])))
    A_p1 = (4π *ρ_cc_p1* r_cc_p1^2)^2 * Λ_cc_p1 * α_w * sqrt(ω_cc_p1)
    F_p1 = (A_p1 / sm.props.dm[k+1]) * (exp(get_p1_dual(sm.props.gamma_turb[k+1])) - exp(get_00_dual(sm.props.gamma_turb[k])))

    # Calculation of all terms for residual 
    mixing_term = (F_p1 - F_00) / dm_face_p1
    omega_var_term = ω_face_00 * dgammadt_face_00 #ω_face_00 *
    source_term = α₁_face_00 * SA_face_00 * sqrt(ω_face_00)
    turb_dissipation_term = C_d * (ω_face_00)^(3/2) / Λ_face_00
    rad_dissipation_term = ω_face_00 / τᵣ_face_00
    excess_term = C_d * (c_s_face_00 * 1e-7)^3 / Λ_face_00  

    return omega_var_term - source_term + turb_dissipation_term + rad_dissipation_term - excess_term - mixing_term
end 
end


# function gammaTurb(sm::StellarModel, k::Int)

# ###Constants###
# C_d = 8/3 * sqrt(2/3)
# α_w = 0.25 
#   ###
# if k == 1
#     γ_face_00 = get_00_dual(sm.props.gamma_turb[k])
#     ω_face_00 = exp(γ_face_00)
#     #dm_cell_00 = sm.props.dm[k] 
#     m_cell_00 = sm.props.m[k]
#     r_face_00 = exp(get_00_dual(sm.props.lnr[k]))
#     P_face_00 = exp(get_00_dual(sm.props.lnP_face[k]))
#     ρ_face_00 = exp(get_00_dual(sm.props.lnρ_face[k]))
#     T_face_00 = exp(get_00_dual(sm.props.lnT_face[k]))
#     κ_face_00 = get_00_dual(sm.props.κ_face[k])
#     ∇ₐ_face_00 = get_00_dual(sm.props.∇ₐ_face[k])
#     cₚ_face_00 = get_00_dual(sm.props.cₚ_face[k])
#     L_face_00 = get_00_dual(sm.props.L[k]) * LSUN

#     Hₚ_face_00 = P_face_00 / (ρ_face_00 * m_cell_00 * CGRAV / r_face_00^2)
#     Λ_face_00 = 1 / (1 / Hₚ_face_00 + 1 / r_face_00)
#     ∇ᵣ_face_00 = (3 * κ_face_00 * L_face_00 * P_face_00) / (16π * CRAD * CLIGHT * CGRAV * m_cell_00 * T_face_00^4)
#     τᵣ_face_00 = (cₚ_face_00 * κ_face_00 * ρ_face_00^2 * Λ_face_00^2) / (48 * SIGMA_SB * T_face_00^3)
#     c_s_face_00 = sqrt(P_face_00 / ρ_face_00)
#     k_rad_face_00 = (16 * SIGMA_SB * T_face_00^3) / (3 * κ_face_00 * ρ_face_00)  ##
#     α₂_face_00 = ρ_face_00 * cₚ_face_00 * 0.5 * sqrt(2 / 3) * Λ_face_00 * sqrt(ω_face_00)
#     α₁_face_00 = ∇ₐ_face_00 * T_face_00 * Λ_face_00 * 0.5 * sqrt(2 / 3) * cₚ_face_00 / Hₚ_face_00^2
#     SA_face_00 = (∇ᵣ_face_00 - ∇ₐ_face_00) * (1 + α₂_face_00 / k_rad_face_00)^(-1)
#     dgammadt_face_00 = (γ_face_00- get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt
#     # Different terms for residual at k = 1

#     # mixing_term =  (F_p1/ dm_cell_00)
#     omega_var_term = ω_face_00 * dgammadt_face_00  #ω_face_00 *
#     source_term = α₁_face_00 * SA_face_00 * sqrt(ω_face_00)
#     turb_dissipation_term = C_d * (ω_face_00)^(3/2) / Λ_face_00
#     rad_dissipation_term = ω_face_00 / τᵣ_face_00
#     excess_term = C_d * (c_s_face_00 * 1e-7)^3 / Λ_face_00

#     return omega_var_term - source_term + turb_dissipation_term + rad_dissipation_term - excess_term #- mixing_term
# end

# if k == sm.props.nz

#     γ_cc_00 = get_00_dual(sm.props.gamma_turb[k]) 
#     ω_cc_00 = exp(γ_cc_00) ###
#     #dm_face_00 = 0.5*(sm.props.dm[k] + sm.props.dm[k-1])
#     m_face_00 = sm.props.m[k]
#     r_cc_00 = exp(get_00_dual(sm.props.lnr[k]))
#     P_cc_00 = get_00_dual(sm.props.eos_res[k].P)
#     ρ_cc_00 = get_00_dual(sm.props.eos_res[k].ρ)
#     T_cc_00 = get_00_dual(sm.props.eos_res[k].T)
#     κ_cc_00 = get_00_dual(sm.props.κ[k])
#     ∇ₐ_cc_00 = get_00_dual(sm.props.eos_res[k].∇ₐ)
#     cₚ_cc_00 = get_00_dual(sm.props.eos_res[k].cₚ)
#     L_cc_00 = get_00_dual(sm.props.L[k])* LSUN

#     Hₚ_cc_00 = P_cc_00 / (ρ_cc_00 * m_face_00 * CGRAV / r_cc_00^2)
#     Λ_cc_00 = 1 / (1 / Hₚ_cc_00 + 1 / r_cc_00)
#     ∇ᵣ_cc_00 = get_00_dual(sm.props.turb_res[k].∇ᵣ)
#     τᵣ_cc_00 = (cₚ_cc_00 * κ_cc_00 * ρ_cc_00^2 * Λ_cc_00^2) / (48 * SIGMA_SB * T_cc_00^3)
#     c_s_cc_00 = sqrt(P_cc_00 / ρ_cc_00)
#     k_rad_cc_00 = (16 * SIGMA_SB * T_cc_00^3) / (3 * κ_cc_00 * ρ_cc_00)  ##
#     α₂_cc_00 = ρ_cc_00 * cₚ_cc_00 * 0.5 * sqrt(2 / 3) * Λ_cc_00 * sqrt(ω_cc_00)
#     α₁_cc_00 = ∇ₐ_cc_00 * T_cc_00 * Λ_cc_00 * 0.5 * sqrt(2 / 3) * cₚ_cc_00 / Hₚ_cc_00^2
#     SA_cc_00 = (∇ᵣ_cc_00 - ∇ₐ_cc_00) * (1 + α₂_cc_00 / k_rad_cc_00)^(-1)    
#     dgammadt_cc_00 = (γ_cc_00- (get_value(sm.start_step_props.gamma_turb[k]) )) / sm.props.dt


#     omega_var_term = ω_cc_00* dgammadt_cc_00 #ω_cc_00*
#     source_term = α₁_cc_00 * SA_cc_00 * sqrt(ω_cc_00)
#     turb_dissipation_term = C_d * (ω_cc_00)^(3/2) / Λ_cc_00
#     rad_dissipation_term = ω_cc_00 / τᵣ_cc_00
#     excess_term = C_d * (c_s_cc_00 * 1e-7)^3 / Λ_cc_00   
#     return omega_var_term - source_term + turb_dissipation_term + rad_dissipation_term - excess_term  #- mixing_term


# end 

# begin
#     γ_face_00 = get_00_dual(sm.props.gamma_turb[k])
#     ω_face_00 = exp(γ_face_00)  ###
#     #dm_cell_00 = sm.props.dm[k] 
#     m_cell_00 = sm.props.m[k]
#     r_face_00 = exp(get_00_dual(sm.props.lnr[k]))
#     P_face_00 = exp(get_00_dual(sm.props.lnP_face[k]))
#     ρ_face_00 = exp(get_00_dual(sm.props.lnρ_face[k]))
#     T_face_00 = exp(get_00_dual(sm.props.lnT_face[k]))
#     κ_face_00 = (get_00_dual(sm.props.κ_face[k]))
#     ∇ₐ_face_00 = get_00_dual(sm.props.∇ₐ_face[k])
#     cₚ_face_00 = get_00_dual(sm.props.cₚ_face[k])
#     L_face_00 = get_00_dual(sm.props.L[k]) * LSUN

#     Hₚ_face_00 = P_face_00 / (ρ_face_00 * m_cell_00 * CGRAV / r_face_00^2)
#     Λ_face_00 = 1 / (1 / Hₚ_face_00 + 1 / r_face_00)
#     ∇ᵣ_face_00 = (3 * κ_face_00 * L_face_00 * P_face_00) / (16π * CRAD * CLIGHT * CGRAV * m_cell_00 * T_face_00^4)
#     τᵣ_face_00 = (cₚ_face_00 * κ_face_00 * ρ_face_00^2 * Λ_face_00^2) / (48 * SIGMA_SB * T_face_00^3)
#     c_s_face_00 = sqrt(P_face_00 / ρ_face_00)
#     k_rad_face_00 = (16 * SIGMA_SB * T_face_00^3) / (3 * κ_face_00 * ρ_face_00)  ##
#     α₂_face_00 = ρ_face_00 * cₚ_face_00 * 0.5 * sqrt(2 / 3) * Λ_face_00 * sqrt(ω_face_00)
#     α₁_face_00 = ∇ₐ_face_00 * T_face_00 * Λ_face_00 * 0.5 * sqrt(2 / 3) * cₚ_face_00 / Hₚ_face_00^2
#     SA_face_00 = (∇ᵣ_face_00 - ∇ₐ_face_00) * (1 + α₂_face_00 / k_rad_face_00)^(-1)
#     dgammadt_face_00 = (γ_face_00- get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt
#     # Different terms for residual at k = 1

#     # mixing_term =  (F_p1/ dm_cell_00)
#     omega_var_term = ω_face_00 * dgammadt_face_00  #ω_face_00 *
#     source_term = α₁_face_00 * SA_face_00 * sqrt(ω_face_00)
#     turb_dissipation_term = C_d * (ω_face_00)^(3/2) / Λ_face_00
#     rad_dissipation_term = ω_face_00 / τᵣ_face_00
#     excess_term = C_d * (c_s_face_00 * 1e-7)^3 / Λ_face_00

#     return omega_var_term - source_term + turb_dissipation_term + rad_dissipation_term - excess_term #- mixing_term
# end

# end


#  function gammaTurb(sm::StellarModel, k::Int)
    
    
    
#     γ₀ = get_00_dual(sm.props.gamma_turb[k])
#     ω = exp(γ₀)
#     r₀ = exp(get_00_dual(sm.props.lnr[k]))
#     m₀ = sm.props.m[k]
#     ∇ₐ = get_00_dual(sm.props.∇ₐ_face[k])
#     cₚ =  get_00_dual(sm.props.cₚ_face[k])
#     κ = get_00_dual(sm.props.κ_face[k])

#     ρ_face = exp(get_00_dual(sm.props.lnρ_face[k]))
#     P_face = exp(get_00_dual(sm.props.lnP_face[k]))
#     T_face = exp(get_00_dual(sm.props.lnT_face[k]))



#      if k == sm.props.nz
#         ρ₀ = get_00_dual(sm.props.eos_res[k].ρ)
#         P₀ = get_00_dual(sm.props.eos_res[k].P)
#         r₀ = exp(get_00_dual(sm.props.lnr[k]))
#         ∇ₐ = get_00_dual(sm.props.eos_res[k].∇ₐ)
#         ∇ᵣ = get_00_dual(sm.props.turb_res[k].∇ᵣ)
#         cₚ = get_00_dual(sm.props.eos_res[k].cₚ)
#         κ = get_00_dual(sm.props.κ[k])
#         Hₚ = P₀ / (ρ₀ * CGRAV * m₀ / r₀^2)
#         T₀ = get_00_dual(sm.props.eos_res[k].T)
#         Λ = 1/(1/Hₚ + 1/r₀) 
#         τᵣ = cₚ * κ * ρ₀^2 * Λ^2 / (48 * SIGMA_SB * T₀^3 )
#         c_s  = sqrt(P₀/ρ₀)
#         k_rad = 16 * SIGMA_SB * T₀^3 / (3 * κ * ρ₀)
#         α₁ = ∇ₐ * T₀ * Λ * 0.5*sqrt(2/3) * cₚ/ Hₚ^2
#         C_d = 8/3 * sqrt(2/3)
#         α₂ = ρ₀*cₚ*0.5*sqrt(2/3)*Λ*sqrt(ω)
#         SA = (∇ᵣ - ∇ₐ)*(1 + α₂/k_rad)^(-1)
#         ∇ = ∇ₐ + SA
#         dgammadt = (γ₀ - get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt
#         @show  cₚ.value, ∇ₐ.value, ∇ᵣ.value,τᵣ.value
#         return ω*dgammadt - α₁* SA * sqrt(ω) + C_d *ω^(3/2)/Λ + ω/τᵣ - C_d*(c_s*1e-7)^3/Λ
#     end

#     # for k != nz 
#     L = get_00_dual(sm.props.L[k]) * LSUN
#     Hₚ = P_face / (ρ_face * CGRAV * m₀ / r₀^2) #defined at face 
#     Λ = 1/(1/Hₚ + 1/r₀) 
#     τᵣ = cₚ * κ * ρ_face^2 * Λ^2 / (48 * SIGMA_SB * T_face^3 )
#     c_s  = sqrt(P_face/ρ_face)
#     k_rad = 16 * SIGMA_SB * T_face^3 / (3 * κ * ρ_face)
#     α₁ = ∇ₐ * T_face * Λ * 0.5*sqrt(2/3) * cₚ/ Hₚ^2
#     C_d = 8/3 * sqrt(2/3)
#     α₂ = ρ_face*cₚ*0.5*sqrt(2/3)*Λ*sqrt(ω)
#     ∇ᵣ = 3 * κ * L * P_face / (16π * CRAD * CLIGHT * CGRAV * m₀ * T_face^4)
#     SA = (∇ᵣ - ∇ₐ)*(1 + α₂/k_rad)^(-1)
#     ∇ = ∇ₐ + SA
    
#     dgammadt = (γ₀ - get_value(sm.start_step_props.gamma_turb[k])) / sm.props.dt
#     # if (k== 1)
#     #     @show (α₁* SA * sqrt(ω)), (ω*dgammadt), (C_d *ω^(3/2)/Hₚ), (ω/τᵣ),  (C_d*(c_s*1e-10)^3/Hₚ)
#     #     @show (ω).value, γ₀.value, ρ_face.value, P_face.value, T_face.value
#     #  end
#     return  ω*dgammadt - α₁* SA * sqrt(ω) + C_d *ω^(3/2)/Λ + ω/τᵣ - C_d*(c_s*1e-7)^3/Λ
# end