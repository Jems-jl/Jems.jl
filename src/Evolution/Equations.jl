using Jems.DualSupport

"""
    equationHSE(sm::StellarModel, k::Int,
                varm1::Matrix{TT}, var00::Matrix{TT}, varp1::Matrix{TT},
                eosm1::EOSResults{TT}, eos00::EOSResults{TT}, eosp1::EOSResults{TT},
                rates::Matrix{TT},
                κm1::TT, κ00::TT, κp1::TT)::TT where {TT<:Real}

Default equation of hydrostatic equilibrium. Evaluates for cell `k` of StellarModel `sm` to what degree hydrostatic
equilibrium is satisfied.

# Arguments

  - `sm`: Stellar Model
  - `k`: cell number to consider
  - `varm1`: Matrix holding the dual numbers of the previous cell (`k-1`)
  - `var00`: Matrix holding the dual numbers of this cell (`k`)
  - `varp1`: Matrix holding the dual numbers of the next cell (`k+1`)
  - `eosm1`: EOSResults object holding the results of the EOS evaluation of the previous cell (`k-1`)
  - `eos00`: EOSResults object holding the results of the EOS evaluation of the current cell (`k`)
  - `eosp1`: EOSResults object holding the results of the EOS evaluation of the next cell (`k+1`)
  - `κm1`: Opacity evaluated at the previous cell (`k-1`)
  - `κ00`: Opacity evaluated at the current cell (`k`)
  - `κp1`: Opacity evaluated at the next cell (`k+1`)

# Returns

Residual of comparing dlnP/dm with -GM/4πr^4, where the latter is evaluated at the face of cell `k` and `k+1`.
"""
function equationHSE(sm::StellarModel, k::Int)
    if k == sm.nz  # atmosphere boundary condition
        lnP₀ = get_00_dual(sm.props.eos_res[k].lnP)
        r₀ = exp(get_00_dual(sm.props.lnr[k]))
        g₀ = CGRAV * sm.mstar / r₀^2
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

    lnPface = (sm.dm[k+1] * lnP₀ + sm.dm[k] * lnP₊) / (sm.dm[k] + sm.dm[k + 1])
    r₀ = exp(get_00_dual(sm.props.lnr[k]))
    dm = 0.5*(sm.dm[k + 1] + sm.dm[k])

    return (exp(lnPface) * (lnP₊ - lnP₀) / dm + CGRAV * sm.m[k] / (4π * r₀^4)) /
           (CGRAV * sm.m[k] / (4π * r₀^4))
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

Residual of comparing dlnT/dm with -∇*GMT/4πr^4P, where the latter is evaluation at the face of cell `k` and `k+1`.
"""
function equationT(sm::StellarModel, k::Int)
    lnT₀ = get_00_dual(sm.props.eos_res[k].lnT)
    if k == sm.nz  # atmosphere boundary condition
        L₀ = get_00_dual(sm.props.L[k]) * LSUN
        r₀ = exp(get_00_dual(sm.props.lnr[k]))
        return lnT₀ - log(L₀ / (BOLTZ_SIGMA * 4π * r₀^2)) / 4  # Eddington gray, ignoring radiation pressure term
    end
    κ00 = get_00_dual(sm.props.κ[k])
    κp1 = get_p1_dual(sm.props.κ[k+1])
    κface = exp((sm.dm[k] * log(κ00) + sm.dm[k + 1] * log(κp1)) / (sm.dm[k] + sm.dm[k + 1]))
    L₀ = get_00_dual(sm.props.L[k]) * LSUN
    r₀ = exp(get_00_dual(sm.props.lnr[k]))

    lnP₀ = get_00_dual(sm.props.eos_res[k].lnP)
    lnP₊ = get_p1_dual(sm.props.eos_res[k+1].lnP)
    Pface = exp((sm.dm[k] * lnP₀ + sm.dm[k + 1] * lnP₊) /
                (sm.dm[k] + sm.dm[k + 1]))
    lnT₊ = get_p1_dual(sm.props.eos_res[k+1].lnT)
    Tface = exp((sm.dm[k] * lnT₀ + sm.dm[k + 1] * lnT₊) / (sm.dm[k] + sm.dm[k + 1]))

    ∇ᵣ = 3κface * L₀ * Pface / (16π * CRAD * CLIGHT * CGRAV * sm.m[k] * Tface^4)
    ∇ₐ_p1 = get_p1_dual(sm.props.eos_res[k+1].∇ₐ)
    ∇ₐ_00 = get_00_dual(sm.props.eos_res[k].∇ₐ)
    ∇ₐ = (sm.dm[k] * ∇ₐ_00 + sm.dm[k + 1] * ∇ₐ_p1) / (sm.dm[k] + sm.dm[k + 1])

    if (∇ᵣ < ∇ₐ)
        return (Tface * (lnT₊ - lnT₀) / sm.dm[k] +
                CGRAV * sm.m[k] * Tface / (4π * r₀^4 * Pface) * ∇ᵣ) /
               (CGRAV * sm.m[k] * Tface / (4π * r₀^4 * Pface))  # only radiative transport
    else  # should do convection here
        return (Tface * (lnT₊ - lnT₀) / sm.dm[k] +
                CGRAV * sm.m[k] * Tface / (4π * r₀^4 * Pface) * ∇ₐ) /
               (CGRAV * sm.m[k] * Tface / (4π * r₀^4 * Pface))  # only radiative transport
    end
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
    dTdt = (T₀ - exp(sm.ssi.lnT[k])) / sm.ssi.dt
    dPdt = (P₀ - exp(sm.ssi.lnP[k])) / sm.ssi.dt

    ϵnuc::typeof(L₀) = 0
    for i in eachindex(sm.network.reactions)
        ϵnuc += get_00_dual(sm.props.rates[k,i])*sm.network.reactions[i].Qvalue
    end
    if k > 1
        L₋ = get_m1_dual(sm.props.L[k-1]) * LSUN
        return ((L₀ - L₋) / sm.dm[k] - ϵnuc + cₚ * dTdt - (δ / ρ₀) * dPdt)  # no neutrinos
    else
        return (L₀ / sm.dm[k] - ϵnuc + cₚ * dTdt - (δ / ρ₀) * dPdt)  # no neutrinos
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
        actual_dr³_dm = (r₀^3 - r₋^3) / sm.dm[k]
    else
        actual_dr³_dm = (r₀^3) / sm.dm[k]
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

TBD

# Returns

Residual of comparing dX_i/dt with its computed reaction rate
"""
function equation_composition(sm::StellarModel, k::Int, iso_name::Symbol)
    # Get mass fraction for this iso
    X = get_00_dual(sm.props.xa[k, sm.network.xa_index[iso_name]])

    dXdt_nuc::typeof(X) = 0
    reactions_in = sm.network.species_reactions_in[sm.network.xa_index[iso_name]]
    for reaction_in in reactions_in
        rate = get_00_dual(sm.props.rates[k,reaction_in[1]])
        dXdt_nuc = dXdt_nuc - rate*reaction_in[2]*Chem.isotope_list[iso_name].A*AMU
    end
    reactions_out = sm.network.species_reactions_out[sm.network.xa_index[iso_name]]
    for reaction_out in reactions_out
        rate = get_00_dual(sm.props.rates[k,reaction_out[1]])
        dXdt_nuc = dXdt_nuc + rate*reaction_out[2]*Chem.isotope_list[iso_name].A*AMU
    end

    Xi = sm.ssi.ind_vars[(k - 1) * sm.nvars + sm.vari[iso_name]]

    return (X - Xi) / sm.ssi.dt - dXdt_nuc
end
