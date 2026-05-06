using ForwardDiff
using Interpolations

"""
    calc_structure_residuals(lny_guess, lny_prev, P_prev ,∇_ad_prev, m_face, dm_cell, dm_cell_p1,  dm_face, sm, xa, species_names, is_core)

Evaluates the numerical residuals for the stellar structure equations (Continuity, Hydrostatic Equilibrium, 
and an Adiabatic Temperature Gradient) between two adjacent mass cells.
Returns an array of residuals `[F_r, F_P, F_T]` which should approach zero upon convergence.
"""
function calc_structure_residuals(ln_state_next :: AbstractVector{T}, ln_state_prev, P_prev ,∇_ad_prev, m_face, dm_next, dm_prev,  dm_face, sm, xa, species_names, is_core) where {T <: Real}
    # Next Cell Variables 
    r_next = exp(ln_state_next[1])
    ρ_next = exp(ln_state_next[2])
    T_next = exp(ln_state_next[3])
    
    # Current cell variables calculated in last iteration 
    r_prev = exp(ln_state_prev[1])
    T_prev = exp(ln_state_prev[3])
    ρ_prev = exp(ln_state_prev[2])
    
    # Dual number handling for ForwardDiff
    xa_dual = [T(x) for x in xa] 
    r_eos = EOSResults{T}()

    # EOS for cell k+1 
    set_EOS_resultsTρ!(sm.eos, r_eos, ln_state_next[3], ln_state_next[2], xa_dual, species_names)
    P_next = r_eos.P
    ∇_ad_next = r_eos.∇ₐ 
    # Face average of some values 
    ln_P_face = (dm_prev * log(P_prev) + dm_next * log(P_next)) / (dm_prev + dm_next)
    ∇_ad_face = (dm_prev * ∇_ad_prev + dm_next * ∇_ad_next) / (dm_prev + dm_next) 

    # Evaluate Residuals
    # --- Continuity Equation ---  
    # Expected dr^3 / dm based on the guessed density
    expected_dr3_dm = 3.0 / (4.0 * π * ρ_next)
    actual_dr3_dm = (r_next^3 - r_prev^3) / dm_next
    F_r = (expected_dr3_dm - actual_dr3_dm) * ρ_next # Scaled for numerical stability
    
    # --- Hydrostatic Equilibrium & Temperature Gradient ---
    if is_core
        # Core expansion: g goes to 0, use volume-averaged rmid
        rmid = 0.5 * (r_prev + r_next)
        F_P = 1.0 - (P_next + CGRAV * ρ_prev^2 * (2.0 * π / 3.0) * rmid^2) / P_prev
        
        dlnP = log(P_next) - log(P_prev)
        dlnT = ln_state_next[3] - ln_state_prev[3]
        F_T = dlnT - dlnP * ∇_ad_face
    else
        # Standard envelope finite difference
        dlnP = log(P_next) - log(P_prev)
        dlnT = ln_state_next[3] - ln_state_prev[3] 
        
        dlnP_expected = - (CGRAV * m_face * dm_face) / (4.0 * π * r_prev^4 * exp(ln_P_face))
        dlnT_expected = dlnP_expected * ∇_ad_face
        
        F_P = dlnP - dlnP_expected
        F_T = dlnT - dlnT_expected
    end
    
    return [F_r, F_P, F_T]
end

"""
    step_integrator(lny_prev, P_prev, ∇_ad_prev, m_face_val, dm_cell, dm_cell_p1, dm_face_val, sm, xa, species_names, is_core)

Performs a single spatial integration step outward to calculate the state (radius, density, temperature) 
of the next mass cell. Uses an initial heuristic guess followed by a Newton-Raphson root-finding 
algorithm.
"""
function step_integrator(lny_prev, P_prev, ∇_ad_prev, m_face_val, dm_cell,dm_cell_p1,dm_face_val, sm, xa, species_names, is_core)
    
    r_prev = exp(lny_prev[1])
    ρ_prev = exp(lny_prev[2])
    T_prev = exp(lny_prev[3])
    r_guess = (r_prev^3 + (3 * dm_cell) / (4 * π * ρ_prev))^(1/3)

    if is_core
        ΔP_est = -CGRAV * ρ_prev^2 * (2π / 3) * r_prev^2
    else
        ΔP_est = -CGRAV * m_face_val * dm_face_val / (4 * π * r_prev^4)
    end
    P_guess = max(P_prev + ΔP_est, P_prev * 0.5)


    dlnP_est = log(P_guess) - log(P_prev)
    dlnT_est = dlnP_est * (∇_ad_prev)
    T_guess = T_prev * exp(dlnT_est)

    Prad_prev = (CRAD / 3.0) * T_prev^4
    Prad_guess = (CRAD / 3.0) * T_guess^4

    Pgas_prev = max(P_prev - Prad_prev, P_prev * 0.01)
    Pgas_guess = max(P_guess - Prad_guess, P_guess * 0.01)

    ρ_guess = ρ_prev * (Pgas_guess / Pgas_prev) * (T_prev / T_guess)
    ρ_guess = max(ρ_guess, ρ_prev * 0.1)

    lny_guess = [log(r_guess), log(ρ_guess), log(T_guess)]
    
    correction_tolerance = 1e-8
    residual_tolerance = 1e-6

    for i in 1:100
    
        F = calc_structure_residuals(lny_guess, lny_prev, P_prev, ∇_ad_prev, m_face_val, dm_cell, dm_cell_p1, dm_face_val, sm, xa, species_names, is_core)
        
        J = ForwardDiff.jacobian(y -> calc_structure_residuals(y, lny_prev, P_prev, ∇_ad_prev, m_face_val, dm_cell, dm_cell_p1, dm_face_val, sm, xa, species_names, is_core), lny_guess)
        Δy = J \ F
        max_step = maximum(abs.(Δy))
        min_step = minimum(abs.(Δy))
        max_allowed_step = 0.05
        
        if max_step > max_allowed_step 
            Δy = Δy .* (max_allowed_step / max_step)
        end 
        lny_guess = lny_guess - Δy
       
        if min_step < correction_tolerance && maximum(abs.(F)) < residual_tolerance
            return lny_guess
        
        end     
    end
    return lny_guess 
end

"""
    shoot_star(lnρc, lnTc, sm, xa, species_names, nz, dms, m_face)

Integrates the stellar structure equations from the core (`lnρc`, `lnTc`) outwards to the surface 
using the shooting method. Iteratively applies `step_integrator` over the defined mass grid `dms`.
Returns the full radial arrays for log-radius, log-density, log-temperature, pressure, and the adiabatic gradient.
"""
function shoot_star(lnρc, lnTc, sm, xa, species_names, nz, dms, m_face)
    ln_r_arr = zeros(nz)
    ln_ρ_arr = zeros(nz)
    ln_T_arr = zeros(nz)
    P_arr = zeros(nz)
    ∇_ad_arr = zeros(nz)

    # Central condition/values 
    r_c = EOSResults{Float64}()
    set_EOS_resultsTρ!(sm.eos, r_c, lnTc, lnρc, xa, species_names)
    Pc = r_c.P
    ∇_ad_c = r_c.∇ₐ 
    ρc = exp(lnρc)

    r1 = (3 * dms[1] / (4 * π * ρc))^(1/3) # Radius of first cell 

    
    ln_r_arr[1] = log(r1) #log radius of first cell in the array 
    ln_ρ_arr[1] = lnρc # first cell value 
    ln_T_arr[1] = lnTc # first cell value 
    P_arr[1] = Pc #first cell value 
    ∇_ad_arr[1] = ∇_ad_c 
    # Integrate outward
    for k in 2:nz
        lny_prev = [ln_r_arr[k-1], ln_ρ_arr[k-1], ln_T_arr[k-1]]
        dm_cell = dms[k]
        dm_cell_p1 = dms[k-1]
        dm_face_val = 0.5 * (dms[k-1] + dms[k])
        m_face_val = m_face[k-1]

        is_core = (k == 2)
        
        lny_next = step_integrator(lny_prev, P_arr[k-1], ∇_ad_arr[k-1], m_face_val, dm_cell, dm_cell_p1, dm_face_val, sm, xa, species_names, is_core)
        
        ln_r_arr[k], ln_ρ_arr[k], ln_T_arr[k] = lny_next[1], lny_next[2], lny_next[3]
        
        r_k = EOSResults{Float64}()
        set_EOS_resultsTρ!(sm.eos, r_k, ln_T_arr[k], ln_ρ_arr[k], xa, species_names)
        P_arr[k] = r_k.P
        ∇_ad_arr[k] = r_k.∇ₐ
            
        
    end 
    
    return ln_r_arr, ln_ρ_arr, ln_T_arr, P_arr, ∇_ad_arr
end

"""
    n_polytrope_central_initial_condition(n, sm, X, Z, Dfraction, abundanceList, M, R)

Calculates an initial estimate for the central density (`lnρc`) and central temperature (`lnTc`) 
by mapping a user-provided Mass and Radius to a fully realized Lane-Emden polytropic structure.
"""
function n_polytrope_central_initial_condition(n, sm::StellarModel, X, Z, Dfraction, abundanceList::AbundanceList, M::Real, R::Real)
    xvals, yvals, zvals = RungeKutta_LaneEmden(n)
    (θ_n, ξ_1, derivative_θ_n) = (linear_interpolation(xvals,yvals), xvals[end],linear_interpolation(xvals,zvals))
    rn = R / ξ_1  # ξ is defined as r/rn, where rn^2=(n+1)Pc/(4π G ρc^2)
    ρc = M / (4π * rn^3 * (-ξ_1^2 * derivative_θ_n(ξ_1)))
    lnρc = log(ρc)
    Pc = 4π * CGRAV * rn^2 * ρc^2 / (n + 1)
    massfractions = get_mass_fractions(abundanceList, sm.network.species_names, X, Z, Dfraction)
    μ = EOS.get_μ_IdealEOS(collect(values(massfractions)), sm.network.species_names)
    lnT_initial = log(Pc * μ / (CGAS * ρc))
    lnTc = getlnT_NewtonRhapson(lnT_initial, lnρc, Pc, massfractions, sm.eos)
    
    return lnρc, lnTc
end 

"""
    ismonotonic(A :: Array, cmp = >)

Helper function that checks if an array strictly follows a monotonic trend based on the comparator `cmp`.
Defaults to checking for strictly decreasing arrays for verifying physical pressure drops.
"""
function ismonotonic(A :: Array, cmp = >)
    current = A[1]
    for i in 2:(length(A)-1)
        newval = A[i]
        cmp(newval,current) && return false
        current = newval
    end 
    return true 
 end 

 """
    surface_boundary(lnρc, lnTc, sm, xa, species_names, nz, dms, m_face)

Matches the central initial conditions to the surface conditions.
First drops `lnρc` until the resulting pressure profile from `shoot_star` is monotonically decreasing.
Then uses a bisection search on `lnρc` to find the exact core density required for the outward-integrated 
surface radiative luminosity to match the Stefan-Boltzmann law.
"""
 function surface_boundary(lnρc, lnTc, sm, xa, species_names, nz, dms, m_face)
    step = 0.5
    max_iter = 10 
    iter = 0 
    local ln_r_arr, ln_ρ_arr, ln_T_arr, P_arr, ∇_ad_arr
    
    # monotonicity check 
    while iter < max_iter
        ln_r_arr, ln_ρ_arr, ln_T_arr, P_arr, ∇_ad_arr = shoot_star(lnρc, lnTc, sm, xa, species_names, nz, dms, m_face)
        if ismonotonic(P_arr)
            break 
        else 
            lnρc -= step     
        end 
        iter += 1
    end 

    # Bisection solver for surface luminosity
    ρ_low = lnρc  
    ρ_high = lnρc + step 
    residual_threshold = 1e-10

    iter_bisect = 0 
    max_bisect_iter = 50 

    while iter_bisect < max_bisect_iter
        ρ_curr = 0.5 * (ρ_low + ρ_high) 
        ln_r_arr, ln_ρ_arr, ln_T_arr, P_arr, ∇_ad_arr = shoot_star(ρ_curr, lnTc, sm, xa, species_names, nz, dms, m_face)

        if !ismonotonic(P_arr)
            ρ_high = ρ_curr 
            iter_bisect += 1
            continue 
        end

        κ_surf = get_opacity_resultsTρ(sm.opacity,ln_T_arr[end], ln_ρ_arr[end], xa, species_names)
        m_surf, T_surf, P_surf, ∇_ad_surf, r_surf = m_face[end], exp(ln_T_arr[end]), P_arr[end], ∇_ad_arr[end], exp(ln_r_arr[end])
        L_rad = (16 * π * CRAD * CLIGHT * CGRAV * m_surf * T_surf^4 * ∇_ad_surf)/ (3 * κ_surf * P_surf)
        L_sb = 4 * π * SIGMA_SB * r_surf^2 * T_surf^4
        err_L = (L_rad - L_sb)/L_sb


        if abs(err_L) < residual_threshold
            break 
        end 
        
        if err_L > 0 
            ρ_high = ρ_curr
        else 
            ρ_low = ρ_curr 
        end 

        iter_bisect += 1
    end 
    
    return ln_r_arr, ln_ρ_arr, ln_T_arr, P_arr, ∇_ad_arr

end

"""
    adiabatic_initial_condition!(n, sm, nz, X, Z, Dfraction, abundanceList, M, R; initial_dt)

Constructs an initial adiabatic stellar model representing a fully convective pre-main sequence star.
Generates the mass grid, guesses core properties from a polytrope approximation, runs the shooting method 
to satisfy surface boundary conditions, and populates the `StellarModel` state variables for evolution.
"""
function adiabatic_initial_condition!(n, sm::StellarModel, nz::Int, X, Z, Dfraction, abundanceList:: AbundanceList, M::Real, R::Real; initial_dt = 100 * SECYEAR)

    # Setup Grid
    logdqs = zeros(length(sm.props.dm))
    for i in 1:nz
        logdqs[i] = get_logdq(i, nz, -12.0, 0.0, -6.0, 200)
    end
    dqs = 10 .^ logdqs
    dqs[nz+1:end] .= 0  
    dqs = dqs ./ sum(dqs)
    dms = dqs .* M
    m_face = cumsum(dms)
    
    # correct m_center
    m_cell = cumsum(dms)
    for i = 1:nz
        if i == 1
            m_cell[i] = 0
        elseif i != nz
            m_cell[i] = m_cell[i] - 0.5 * dms[i]
        end
    end

    massfractions = get_mass_fractions(abundanceList, sm.network.species_names, X, Z, Dfraction)
    species_names = collect(Symbol, keys(massfractions))
    xa = collect(Float64, values(massfractions))
    lnRho_c_guess,  lnT_c_guess = n_polytrope_central_initial_condition(n, sm, X, Z, Dfraction, abundanceList, M, R)
    final_lnr, final_lnρ, final_lnT, final_P, ∇_ad_profile = surface_boundary(lnRho_c_guess, lnT_c_guess, sm, xa, species_names, nz, dms, m_face)
    # Mapping to the StellarModel independent variables (sm.props.ind_vars)
    sm.props.mstar = M
    sm.props.dm = dms
    sm.props.m = m_face

    for i = 1:nz 
        sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnr]] = final_lnr[i]
        sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnρ]] = final_lnρ[i]
        sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnT]] = final_lnT[i]
        κ = get_opacity_resultsTρ(sm.opacity,final_lnT[i], final_lnρ[i], xa, species_names)
        sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[:lum]] = (∇_ad_profile[i]) * (16π * CRAD * CLIGHT * CGRAV * m_face[i] * (exp(final_lnT[i]))^4)/ (3 * κ * final_P[i] * LSUN)
        for (isotope, massfraction) in massfractions
            sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[isotope]] = massfraction 
        end
    end

    sm.props.time = 0.0
    sm.props.dt = initial_dt
    sm.props.dt_next = initial_dt
    sm.props.model_number = 0
    sm.props.nz = nz
end 