using ForwardDiff
using Interpolations

"""
    get_logdq(k::Int, nz::Int, logdq_low::TT, logdq_high::TT, numregion::Int)::TT where {TT<:Real}

Computes the logarithm mass chunk `logdq` for zone `k` of a profile with total zones `nz`, while keeping in mind to
better resolve the first and last `numregion` zones of the profile. It linearly interpolates the value from the inputs
`logdq_low` and `logdq_high` in these regions, while keeping `logdq_high` in the middle zones.
"""
function get_logdq(k::Int, nz::Int, logdq_center::TT, logdq_mid::TT, logdq_surf::TT, numregion::Int)::TT where {TT<:Real}
    if k <= numregion
        return logdq_center + (k - 1) * (logdq_mid - logdq_center) / numregion
    elseif k < nz - numregion
        return logdq_mid
    else
        return logdq_mid + (logdq_surf - logdq_mid) * (k - (nz - numregion)) / numregion
    end
end
"""
    calc_structure_residuals(lny_guess, lny_prev, P_prev ,∇_ad_prev, m_face, dm_cell, dm_cell_p1,  dm_face, sm, xa, species_names, is_core)

Computes the residual for the Continuity equation, Hydrostatic equation, temperature gradient equation (considering a fully adiabatic
temperature gradient). 'ln_state_prev' contains the radius, density, and temperature of the preceding cell, while 'ln_state_guess' contains the 
values for the current cell. 
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

function RungeKutta_LaneEmden(n)
    # defining the Lane-Emden equation
    dydx(x,y,z,n) = z
    dzdx(x,y,z,n) = -y^n -2*z/x 
    # defining the approximation for small x
    y_smallx(x,n) = 1 - 1/6*x^2 + n/120*x^4 -n*(8*n-5)/1520*x^6
    z_smallx(x,n) = - 1/3*x + n/30*x^3 -3*n*(8*n-5)/760*x^5;

    function endOfLoop!(xvals, yvals, zvals, endIndex)
        slope = (yvals[endIndex-1] - yvals[endIndex-2]) / (xvals[endIndex-1] - xvals[endIndex-2])
        xlast = xvals[endIndex-1] - yvals[endIndex-1] / slope
        xvals[1:endIndex-1] = xvals[1:endIndex-1]; 
        # add last entry
        xvals[endIndex] = xlast
        yvals[endIndex] = 0.0
        zvals[endIndex] = zvals[endIndex-1]
        # manually adding the core boundary conditions
        pushfirst!(yvals,1.0)
        pushfirst!(xvals,0.0)
        pushfirst!(zvals,0.0)
        return (xvals[1:endIndex+1],yvals[1:endIndex+1],zvals[1:endIndex+1])
    end

    Δx = 1e-5
    Δx_min = 1e-11
    nsteps = 10_000_000  # maximum number of steps
    xvals = LinRange(Δx,nsteps*Δx,nsteps); xvals = collect(xvals)  # making xvals a mutable array
    yvals = zeros(nsteps); zvals = zeros(nsteps)
    yvals[1] = y_smallx(Δx,n); zvals[1] = z_smallx(Δx,n)

    i = 2
    while i<=nsteps
        try
            x = xvals[i-1]; y = yvals[i-1]; z = zvals[i-1]
            k₁ = Δx*dydx(x,y,z,n); l₁ = Δx*dzdx(x,y,z,n)
            ynew = y + k₁/2
            if ynew < 0.0
                throw(ErrorException("negative value"))
            end
            k₂ = Δx*dydx(x+Δx/2,ynew,z+l₁/2,n); l₂ = Δx*dzdx(x+Δx/2,ynew,z+l₁/2,n)
            ynew = y+k₂/2
            if ynew < 0.0
                throw(ErrorException("negative value"))
            end
            k₃ = Δx*dydx(x+Δx/2,ynew,z+l₂/2,n); l₃ = Δx*dzdx(x+Δx/2,ynew,z+l₂/2,n)
            ynew = y+k₃
            if ynew < 0.0
                throw(ErrorException("negative value"))
            end
            k₄ = Δx*dydx(x+Δx,ynew,z+l₃,n);l₄ = Δx*dzdx(x+Δx,ynew,z+l₃,n)
            ynew = y+k₁/6+k₂/3+k₃/3+k₄/6
            if ynew < 0.0
                throw(ErrorException("negative value"))
            end
            yvals[i] = ynew  # new y value
            zvals[i] = z+l₁/6+l₂/3+l₃/3+l₄/6  # new z value
            i = i+1
        catch e
            if isa(e, ErrorException)
                Δx = Δx/2
                xvals[i] = xvals[i-1] + Δx  # next xvalue is now at a smaller distance from the previous one
                if Δx < Δx_min
                    xvals, yvals, zvals = endOfLoop!(xvals,yvals,zvals,i)
                    break
                end
            else
                throw(e)
            end
        end
    end
    if i>nsteps
        throw(ArgumentError("not able to converge to first zero of Lane-Emden equation"))
    end
    return xvals, yvals, zvals
end

function getlnT_NewtonRhapson(lnT_initial, lnρ, P, massfractions, eos)
    (species_names, xa) = (collect(Symbol,keys(massfractions)), collect(Float64,values(massfractions)))
    ΔlnPmin = 1e-4
    lnT = lnT_initial
    lnT_dual = ForwardDiff.Dual(lnT_initial,1.0)
    lnρ_dual = ForwardDiff.Dual(lnρ,0.0)
    xa_dual = [ForwardDiff.Dual(xa[i],0.0) for i in eachindex(xa)]
    r = EOSResults{typeof(lnT_dual)}()
    set_EOS_resultsTρ!(eos,r,lnT_dual,lnρ_dual,xa_dual,species_names)
    lnP = log(r.P)
    dlnPdlnT = lnP.partials[1]
    i = 0
    # performing the Newton-Rhapson algorithm
    while abs(log(P) - lnP.value) > ΔlnPmin
        lnT = lnT + (log(P) - lnP.value) / dlnPdlnT  # go to the next guess
        lnT_dual = ForwardDiff.Dual(lnT,1.0)  # setting new lnT_dual
        set_EOS_resultsTρ!(eos,r,lnT_dual,lnρ_dual,xa_dual,species_names)
        lnP = log(r.P)
        dlnPdlnT = lnP.partials[1]
        i = i+1
        if i>100
            throw(ArgumentError("not able to converge to equation of state temperature"))
        end
    end
    return lnT
end

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
 
function ismonotonic(A :: Array, cmp = >)
    current = A[1]
    for i in 2:(length(A)-1)
        newval = A[i]
        cmp(newval,current) && return false
        current = newval
    end 
    return true 
 end 

 
 function surface_boundary(lnρc, lnTc, sm, xa, species_names, nz, dms, m_face)
    step = 0.5
    max_iter = 10 
    iter = 0 
    local ln_r_arr, ln_ρ_arr, ln_T_arr, P_arr, ∇_ad_arr
    while iter < max_iter
        ln_r_arr, ln_ρ_arr, ln_T_arr, P_arr, ∇_ad_arr = shoot_star(lnρc, lnTc, sm, xa, species_names, nz, dms, m_face)
        if ismonotonic(P_arr)
            break 
        else 
            lnρc -= step     
        end 
        iter += 1
    end 

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

        κ_surf = Jems.Opacity.get_opacity_resultsTρ(sm.opacity,ln_T_arr[end], ln_ρ_arr[end], xa, species_names)
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

function newton_cubic_solver(A::T, B::T, C::T, D::T, E::T; max_iter=50, tol_corr=1e-10, tol_res=1e-10) where {T <: Real}
    ω = 1e5
    scale = max(abs(D), 1.0)
    
    for i in 1:max_iter
        residual = A*ω^4 + B*ω^3 + C*ω^2 + D*ω + E
        f_prime = 4*A*ω^3 + 3*B*ω^2 + 2*C*ω + D
        
        if abs(f_prime) < floatmin(T)
            return ω 
        end
        
        # 3. Calculate the Correction: Δω = f(ω) / f'(ω)
        corr = residual / f_prime
        
        # Apply the correction to find the next guess
        ω_new = ω - corr
        
        corr_is_small = abs(corr) < tol_corr * max(abs(ω_new), 1.0)
        
        res_is_small = abs(residual) < tol_res * scale
        
        if corr_is_small && res_is_small
            return ω_new
        end

        ω = ω_new
    end
    
    # If the loop exhausts max_iter without passing the checks, 
    # it returns the last calculated value. 
    return ω
end
function tdc_initial_condition!(n, sm::StellarModel, nz::Int, X, Z, Dfraction, abundanceList:: AbundanceList, M::Real, R::Real; initial_dt = 100 * SECYEAR)

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
        for (isotope, massfraction) in massfractions
                sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[isotope]] = massfraction 
            end
    end

    for i = 1:nz 
        r_eos = EOSResults{Float64}()
        set_EOS_resultsTρ!(sm.eos, r_eos,final_lnT[i], final_lnρ[i], xa, species_names)
        κ = Jems.Opacity.get_opacity_resultsTρ(sm.opacity,final_lnT[i], final_lnρ[i], xa, species_names)      
        H_p = final_P[i]/ (exp(final_lnρ[i])* m_face[i] * CGRAV / (exp(final_lnr[i]))^2)
        Λ = (1/H_p + 1/exp(final_lnr[i]))
        k_rad = (16 * SIGMA_SB * (exp(final_lnT[i]))^3) / (3 * κ * exp(final_lnρ[i]))  
        c_s = sqrt(final_P[i]/exp(final_lnρ[i]))

        if i == 1 
            dlnP = log(final_P[i])
            dlnT = final_lnT[i]
            #sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[:lum]] = (dlnT/dlnP) * (16π * CRAD * CLIGHT * CGRAV * m_face[i] * (exp(final_lnT[i]))^4)/ (3 * κ * final_P[i] * LSUN)
            sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[:lum]] = (∇_ad_profile[i] - 0.01) * (16π * CRAD * CLIGHT * CGRAV * m_face[i] * (exp(final_lnT[i]))^4)/ (3 * κ * final_P[i] * LSUN)
            #a1 = ((∇_ad_profile[i] * exp(final_lnT[i]) * Λ * 0.5 * sqrt(2 / 3) * r_eos.cₚ)/H_p^2)*((dlnT/dlnP)-∇_ad_profile[i])
            a1 = ((∇_ad_profile[i] * exp(final_lnT[i]) * Λ * 0.5 * sqrt(2 / 3) * r_eos.cₚ)/H_p^2)*(-0.01)
            a2 = (exp(final_lnρ[i]) * r_eos.cₚ *  0.5 * sqrt(2 / 3) * Λ)/ k_rad
            a3 = (8/3 * sqrt(2/3))/Λ
            a4 = 1/((r_eos.cₚ * κ * exp(final_lnρ[i])^2 * Λ^2) / (48 * SIGMA_SB * (exp(final_lnT[i]))^3))
            a5 = (8/3 * sqrt(2/3)) * (c_s * 1e-4)^3 / Λ
            #Cubic equation call   
            A = a2 * a3  
            B = a3 + a2 * a4
            C = a4
            D = -a1 - a2*a5
            E = -a5
            ω_sol = newton_cubic_solver(A,B,C,D,E)
            sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[:gamma_turb]] = -1.0
        else   
            dlnP = log(final_P[i]) - log(final_P[i-1])
            dlnT = final_lnT[i] -final_lnT[i-1]
            #a1 = ((∇_ad_profile[i] * exp(final_lnT[i]) * Λ * 0.5 * sqrt(2 / 3) * r_eos.cₚ)/H_p^2)*( (dlnT/dlnP) - ∇_ad_profile[i])
            a1 = ((∇_ad_profile[i] * exp(final_lnT[i]) * Λ * 0.5 * sqrt(2 / 3) * r_eos.cₚ)/H_p^2)*(-0.01)
            a2 = (exp(final_lnρ[i]) * r_eos.cₚ *  0.5 * sqrt(2 / 3) * Λ)/ k_rad
            a3 = (8/3 * sqrt(2/3))/Λ
            a4 = 1/((r_eos.cₚ * κ * exp(final_lnρ[i])^2 * Λ^2) / (48 * SIGMA_SB * (exp(final_lnT[i]))^3))
            a5 = (8/3 * sqrt(2/3)) * (c_s * 1e-4)^3 / Λ
            #Cubic equation call       
            A = a2 * a3  
            B = a3 + a2 * a4
            C = a4
            D = -a1 - a2*a5
            E = -a5
            ω_sol = newton_cubic_solver(A,B,C,D,E)
            #sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[:lum]] = (dlnT/dlnP) * (16π * CRAD * CLIGHT * CGRAV * m_face[i] * (exp(final_lnT[i]))^4)/ (3 * κ * final_P[i] * LSUN)
            sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[:lum]] = (∇_ad_profile[i] - 0.01) * (16π * CRAD * CLIGHT * CGRAV * m_face[i] * (exp(final_lnT[i]))^4)/ (3 * κ * final_P[i] * LSUN)
            sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[:gamma_turb]] = 10.0 # This needs to figured out still! 
        end 
    end 
    sm.props.time = 0.0
    sm.props.dt = initial_dt
    sm.props.dt_next = initial_dt
    sm.props.model_number = 0
    sm.props.nz = nz
end 