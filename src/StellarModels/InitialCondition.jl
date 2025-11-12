# use an n=1 polytrope as initial condition, should be better than a constant density star
using ForwardDiff
using Roots
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
    RungeKutta_LaneEmden(n)

Computes the solution of the Lane-Emden equation for polytropic index `n` until the first zero by returning `xvals`,
containing ξ values; `yvals`, containing the corresponding function values θ_n; and `zvals`, containing the derivative.
This naming convention for x (=independent variable), y (=corresponding solution values) and z (=corresponding
derivative values) is used throughout this function. The Lane-Emden equation is solved by performing the Runge-Kutta
method of order 4. The stepsize is allowed to decrease as the function reaches the first zero. The function takes care
of the core boundary conditions at ξ=0. The last ξ value (i.e. the first zero of the Lane-Emden solution) is calculated
by linearly extrapolating the last two points of the solution. This first zero ξ_1 can be found as `xvals`[end].
"""
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

"""
    getlnT_NewtonRhapson(lnT_initial, lnρ, P, massfractions, eos)

Computes the temperature lnT starting from a density `lnρ`, a pressure `P`, a composition `xa`,
a species `species` and an equation of state `eos`. The equation of state gives us the pressure, 
given a certain temperature and density. The idea is to match this pressure with the given pressure `P`
by fitting the temperature. Starting from an initial guess `lnT_initial`, 
the Newton-Rhapson method is used to converge to the final temperature in an interative way.
Each iteration, a new lnT is computed according to the Newton-Rhapson formula using the derivative dlnP/dlnT.
Next, based on this new lnT, the equation of state returns a new pressure. The algorithm stops when
the difference between the calculated pressure and the given pressure is smaller than a certain threshold.
The temperature at which this occurs is returned.

"""
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


"""
    n_polytrope_initial_condition(n,sm::StellarModel, M::Real, R::Real; initial_dt=100 * SECYEAR)

Initializes the stellar model properties `sm.props` with a mesh of size `nz` with values corresponding to a polytrope of
index `n`, setting `sm.props.m`, `sm.props.dm` and the independent variables `sm.props.ind_vars`, etc. accordingly. Also
sets the initial timestep to be taken, `initial_dt`. It first calls the solution to the Lane-Emden equation for index
`n` and then sets radii, densities, pressures and luminosities.
"""
function n_polytrope_initial_condition!(n, sm::StellarModel, nz::Int, X, Z, Dfraction, abundanceList::AbundanceList, M::Real, R::Real; initial_dt=100 * SECYEAR)
    xvals, yvals, zvals = RungeKutta_LaneEmden(n)
    (θ_n, ξ_1, derivative_θ_n) = (linear_interpolation(xvals,yvals), xvals[end],linear_interpolation(xvals,zvals))
    
    logdqs = zeros(length(sm.props.dm))
    for i in 1:nz
        logdqs[i] = get_logdq(i, nz, -12.0, 0.0, -6.0, 200)
    end
    dqs = 10 .^ logdqs
    dqs[nz+1:end] .= 0  # extra entries beyond nz have no mass
    dqs = dqs ./ sum(dqs)
    dms = dqs .* M
    m_face = cumsum(dms)
    m_cell = cumsum(dms)
    # correct m_center
    for i = 1:nz
        if i == 1
            m_cell[i] = 0
        elseif i != nz
            m_cell[i] = m_cell[i] - 0.5 * dms[i]
        end
    end

    rn = R / ξ_1  # ξ is defined as r/rn, where rn^2=(n+1)Pc/(4π G ρc^2)
    ρc = M / (4π * rn^3 * (-ξ_1^2 * derivative_θ_n(ξ_1)))
    Pc = 4π * CGRAV * rn^2 * ρc^2 / (n + 1)
    ξ_cell = zeros(nz)
    ξ_face = zeros(nz)
    function mfunc(ξ, m)
        return m - 4π * rn^3 * ρc * (-ξ^2 * derivative_θ_n(ξ))
    end

    for i = 1:(nz)
        if i == 1
            ξ_cell[i] = 0
        elseif i == nz
            ξ_cell[i] = ξ_1
        else
            mfunc_anon = ξ -> mfunc(ξ, m_cell[i])
            ξ_cell[i] = find_zero(mfunc_anon, (0, ξ_1), Bisection())
        end
        if i == nz
            ξ_face[i] = ξ_1
        else
            mfunc_anon = ξ -> mfunc(ξ, m_face[i])
            ξ_face[i] = find_zero(mfunc_anon, (0, ξ_1), Bisection())
        end
    end
    mfunc_anon = ξ -> mfunc(ξ, 0.99999*M)

    # set radii, pressure and temperature, and mass fractions
    massfractions = get_mass_fractions(abundanceList, sm.network.species_names, X, Z, Dfraction)
    μ = EOS.get_μ_IdealEOS(collect(values(massfractions)), sm.network.species_names)
    for i = 1:nz
        sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnr]] = log(rn * ξ_face[i])
        if i > 1
            P = Pc * (θ_n(ξ_cell[i]))^(n + 1)
            ρ = ρc * (θ_n(ξ_cell[i]))^(n)
        else
            P = Pc
            ρ = ρc
        end
        lnT_initial = log(P * μ / (CGAS * ρ))  # ideal gas temperature as intial guess
        # fit the temperature using the equation of state
        #lnT = getlnT_NewtonRhapson(lnT_initial, log(ρ),P,[1.0,0],[:H1,:He4],sm.eos)
        lnT = getlnT_NewtonRhapson(lnT_initial, log(ρ), P, massfractions, sm.eos)

        sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnρ]] = log(ρ)
        sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnT]] = lnT

        #set mass fractions
        for (isotope, massfraction) in massfractions
            sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[isotope]] = massfraction
        end
    end

    # set m and dm
    sm.props.mstar = M
    sm.props.dm = dms
    sm.props.m = m_face

    # set luminosity
    for i = 1:nz - 1
        Pface = Pc * (θ_n(ξ_face[i]))^(n + 1)
        ρface = ρc * (θ_n(ξ_face[i]))^(n)
        Tfaceinit = Pface * μ / (CGAS * ρface)
        lnTface = getlnT_NewtonRhapson(log(Tfaceinit),log(ρface), Pface, massfractions,sm.eos)
        Tface = exp(lnTface)
       
        dlnT = sm.props.ind_vars[(i) * sm.nvars + sm.vari[:lnT]] - sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnT]]
        if i != 1
            dlnP = log(Pc * (θ_n(ξ_cell[i+1]))^(n + 1)) - log(Pc * (θ_n(ξ_cell[i]))^(n + 1))
        else
            dlnP = log(Pc * (θ_n(ξ_cell[i+1]))^(n + 1)) - log(Pc)
        end
        κ = get_opacity_resultsTρ(sm.opacity, lnTface, log(ρface) ,collect(Float64,values(massfractions)), collect(Symbol,keys(massfractions)))

        sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[:lum]] = (dlnT / dlnP) *
                                                          (16π * CRAD * CLIGHT * CGRAV * m_face[i] * Tface^4) /
                                                          (3κ * Pface * LSUN)
    end

    #set initial values for gamma_turb
    for i= 1:nz
        sm.props.ind_vars[(i - 1) * sm.nvars + sm.vari[:gamma_turb]] = 10
    end

    # modify special cases, just copy values at edges
    sm.props.ind_vars[(nz - 1) * sm.nvars + sm.vari[:lnρ]] = sm.props.ind_vars[(nz - 2) * sm.nvars + sm.vari[:lnρ]]
    sm.props.ind_vars[(nz - 1) * sm.nvars + sm.vari[:lnT]] = sm.props.ind_vars[(nz - 2) * sm.nvars + sm.vari[:lnT]]
    # sm.props.ind_vars[(nz - 1) * sm.nvars + sm.vari[:lum]] = sm.props.ind_vars[(nz - 2) * sm.nvars + sm.vari[:lum]]
    sm.props.ind_vars[(nz - 1) * sm.nvars + sm.vari[:lum]] = sm.props.ind_vars[(nz - 3) * sm.nvars + sm.vari[:lum]]
    sm.props.ind_vars[(nz - 2) * sm.nvars + sm.vari[:lum]] = sm.props.ind_vars[(nz - 3) * sm.nvars + sm.vari[:lum]]

    sm.props.time = 0.0
    sm.props.dt = initial_dt
    sm.props.dt_next = initial_dt
    sm.props.model_number = 0
    sm.props.nz = nz
end

