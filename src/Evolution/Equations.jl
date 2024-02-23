using Jems.DualSupport

"""

    equationHSE(sm::StellarModel, k::Int)


Default equation of hydrostatic equilibrium. Evaluates for cell `k` of StellarModel `sm` to what degree hydrostatic
equilibrium is satisfied.

# Arguments

  - `sm`: Stellar Model
  - `k`: cell number to consider

Must be identical to [`equationContinuity`](@ref) for compatibility with StellarModels.TypeStableEquation

# Returns

Residual of comparing the discretized ∂lnP/∂m with -GM/4πr^4, where the latter is evaluated at the face of cell `k` and
`k+1`.
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


    lnP₀ = get_00_dual(sm.props.eos_res[k].lnP)  # log pressure at cell center of cell k
    lnP₊ = get_p1_dual(sm.props.eos_res[k+1].lnP)  # log pressure at cell center of cell k+1
    lnPface = get_00_dual(sm.props.lnP_face[k])  # log pressure at the face between k and k+1

    r₀ = exp(get_00_dual(sm.props.lnr[k]))
    dm = 0.5*(sm.props.dm[k + 1] + sm.props.dm[k])

    return (exp(lnPface) * (lnP₊ - lnP₀) / dm + CGRAV * sm.props.m[k] / (4π * r₀^4)) /
           (CGRAV * sm.props.m[k] / (4π * r₀^4))
end

"""

    equationT(sm::StellarModel, k::Int)

Default equation of energy transport, evaluated for cell `k` of StellarModel `sm`.


# Arguments

  - `sm`: Stellar Model
  - `k`: cell number to consider

Must be identical to [`equationContinuity`](@ref) for compatibility with StellarModels.TypeStableEquation

# Returns

Residual of comparing the discretized ∂lnT/∂m with -∇*GMT/4πr^4P, where the latter is evaluated at the face of cell `k`
and `k+1`.
"""
function equationT(sm::StellarModel, k::Int)
    lnT₀ = get_00_dual(sm.props.eos_res[k].lnT)
    if k == sm.props.nz  # atmosphere boundary condition
        L₀ = get_00_dual(sm.props.L[k]) * LSUN
        r₀ = exp(get_00_dual(sm.props.lnr[k]))
        return lnT₀ - log(L₀ / (BOLTZ_SIGMA * 4π * r₀^2)) / 4  # Eddington gray, ignoring radiation pressure term
    end
    r₀ = exp(get_00_dual(sm.props.lnr[k]))
    lnT₀ = get_00_dual(sm.props.lnT[k])
    lnT₊ = get_p1_dual(sm.props.lnT[k+1])

    Pface = exp(get_00_dual(sm.props.lnP_face[k]))
    Tface = exp(get_00_dual(sm.props.lnT_face[k]))

    ∇ᵣ = get_00_dual(sm.props.∇ᵣ_face[k])
    ∇ₐ = get_00_dual(sm.props.∇ₐ_face[k])

    if (∇ᵣ < ∇ₐ)
        return (Tface * (lnT₊ - lnT₀) / sm.props.dm[k] +
                CGRAV * sm.props.m[k] * Tface / (4π * r₀^4 * Pface) * ∇ᵣ) /
               (CGRAV * sm.props.m[k] * Tface / (4π * r₀^4 * Pface))  # only radiative transport
    else  # should do convection here
        return (Tface * (lnT₊ - lnT₀) / sm.props.dm[k] +
                CGRAV * sm.props.m[k] * Tface / (4π * r₀^4 * Pface) * ∇ₐ) /
               (CGRAV * sm.props.m[k] * Tface / (4π * r₀^4 * Pface))  # only radiative transport
    end
end

"""

    equationLuminosity(sm::StellarModel, k::Int)

Default equation of energy generation, evaluated for cell `k` of StellarModel `sm`.

# Arguments

  - `sm`: Stellar Model
  - `k`: cell number to consider

Must be identical to [`equationContinuity`](@ref) for compatibility with StellarModels.TypeStableEquation

# Returns

Residual of comparing the discretized ∂L/∂m with ϵnuc - cₚ * ∂T/∂t - (δ / ρ) * ∂P/∂t
"""
function equationLuminosity(sm::StellarModel, k::Int)
    L₀ = get_00_dual(sm.props.L[k]) * LSUN
    ρ₀ = get_00_dual(sm.props.eos_res[k].ρ)
    cₚ = get_00_dual(sm.props.eos_res[k].cₚ)
    δ = get_00_dual(sm.props.eos_res[k].δ)
    T₀ = get_00_dual(sm.props.eos_res[k].T)
    P₀ = get_00_dual(sm.props.eos_res[k].P)
    ∂T∂t = (T₀ - get_cell_value(sm.start_step_props.eos_res[k].T)) / sm.props.dt
    ∂P∂t = (P₀ - get_cell_value(sm.start_step_props.eos_res[k].P)) / sm.props.dt

    ϵnuc::typeof(L₀) = 0
    for i in eachindex(sm.network.reactions)
        ϵnuc += get_00_dual(sm.props.rates[k,i])*sm.network.reactions[i].Qvalue
    end
    if k > 1
        L₋ = get_m1_dual(sm.props.L[k-1]) * LSUN
        return ((L₀ - L₋) / sm.props.dm[k] - ϵnuc + cₚ * ∂T∂t - (δ / ρ₀) * ∂P∂t) / (L₀ / sm.props.m[k])  # no neutrinos
    else
        return (L₀ / sm.props.dm[k] - ϵnuc + cₚ * ∂T∂t - (δ / ρ₀) * ∂P∂t) / (L₀ / sm.props.m[k])  # no neutrinos
    end
end

"""

    equationContinuity(sm::StellarModel, k::Int)

Default equation of mass continuity, evaluated for cell `k` of StellarModel `sm`.

# Arguments

  - `sm`: Stellar Model
  - `k`: cell number to consider

Must be identical to other equations for compatibility with StellarModels.TypeStableEquation

# Returns

Residual of comparing the discretized ∂r^3/∂m with 3/(4πρ)
"""
function equationContinuity(sm::StellarModel, k::Int)
    ρ₀ = get_00_dual(sm.props.eos_res[k].ρ)
    r₀ = exp(get_00_dual(sm.props.lnr[k]))
    expected_∂r³∂m = 3 / (4π * ρ₀)

    if k > 1  # get inner radius
        r₋ = exp(get_m1_dual(sm.props.lnr[k-1]))
        actual_∂r³∂m = (r₀^3 - r₋^3) / sm.props.dm[k]
    else
        actual_∂r³∂m = (r₀^3) / sm.props.dm[k]
    end

    return (expected_∂r³∂m - actual_∂r³∂m) * ρ₀  # times ρ to make eq. dim-less
end

"""

    equation_composition(sm::StellarModel, k::Int, iso_name::Symbol)

Default equation for composition evolution for isotope `iso_name`, evaluated for cell `k` of StellarModel `sm`.

# Arguments

  - `sm`: Stellar Model
  - `k`: cell number to consider

Must be identical to [`equationContinuity`](@ref) for compatibility with StellarModels.TypeStableEquation

# Returns

Residual of comparing the discretized ∂X_i/∂t with its computed reaction rate
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

    Xi = get_cell_value(sm.start_step_props.xa[k, sm.network.xa_index[iso_name]])  # is never a dual!!

    return (X - Xi) / sm.props.dt - dXdt_nuc
end
