# use an n=1 polytrope as initial condition, should be better than a constant density star
using ForwardDiff
using Roots

"""
    θ_n(ξ) = sin(ξ) / ξ, the sinc function
"""
θ_n(ξ) = sin(ξ) / ξ

# create a profile for composition that better resolves edges
"""
    get_logdq(k::Int, nz::Int, logdq_low::TT, logdq_high::TT, numregion::Int)::TT where {TT<:Real}

Computes the logarithm mass chunk `logdq` for zone `k` of a profile with total zones `nz`, while keeping in mind to
better resolve the first and last `numregion` zones of the profile. It linearly interpolates the value from the inputs
`logdq_low` and `logdq_high` in these regions, while keeping `logdq_high` in the middle zones.
"""
function get_logdq(k::Int, nz::Int, logdq_low::TT, logdq_high::TT, numregion::Int)::TT where {TT<:Real}
    if k <= numregion
        return logdq_low + (k - 1) * (logdq_high - logdq_low) / numregion
    elseif k < nz - numregion
        return logdq_high
    else
        return logdq_high + (logdq_low - logdq_high) * (k - (nz - numregion)) / numregion
    end
end

"""
    n1_polytrope_initial_condition(sm::StellarModel, M::Real, R::Real; initial_dt=100 * SECYEAR)

Initializes a stellar model `sm` with values corresponding to an n=1 polytrope, setting the independent variables
`sm.ind_vars`, etc. accordingly. Also sets the initial timestep to be taken, `initial_dt`.
"""
function n1_polytrope_initial_condition!(sm::StellarModel, M::Real, R::Real; initial_dt=100 * SECYEAR)
    logdqs = get_logdq.(1:(sm.nz), sm.nz, -3.0, 0.0, 100)
    dqs = 10 .^ logdqs
    dqs = dqs ./ sum(dqs)
    dms = dqs .* M
    m_face = cumsum(dms)
    m_cell = cumsum(dms)
    # correct m_center
    for i = 1:(sm.nz)
        if i == 1
            m_cell[i] = 0
        elseif i != sm.nz
            m_cell[i] = m_cell[i] - 0.5 * dms[i]
        end
    end

    n = 1  # n = 1 polytrope after all...
    rn = R / π  # ξ is defined as r/rn, where rn^2=(n+1)Pc/(4π G ρc^2)

    ρc = M / (4π * rn^3 * (-π^2 * ForwardDiff.derivative(θ_n, π)))
    Pc = 4π * CGRAV * rn^2 * ρc^2 / (n + 1)

    ξ_cell = zeros(sm.nz)
    ξ_face = zeros(sm.nz)
    function mfunc(ξ, m)
        return m - 4π * rn^3 * ρc * (-(-sin(ξ) + cos(ξ) * ξ))
    end

    for i = 1:(sm.nz)
        if i == 1
            ξ_cell[i] = 0
        elseif i == sm.nz
            ξ_cell[i] = π
        else
            mfunc_anon = ξ -> mfunc(ξ, m_cell[i])
            ξ_cell[i] = find_zero(mfunc_anon, (0, π), Bisection())
        end
        if i == sm.nz
            ξ_face[i] = π
        else
            mfunc_anon = ξ -> mfunc(ξ, m_face[i])
            ξ_face[i] = find_zero(mfunc_anon, (0, π), Bisection())
        end
    end

    # set radii, pressure and temperature, assuming ideal gas without Prad
    for i = 1:(sm.nz)
        μ = 0.5
        XH = 1.0
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnr]] = log(rn * ξ_face[i])
        if i > 1
            P = Pc * (θ_n(ξ_cell[i]))^(n + 1)
            ρ = ρc * (θ_n(ξ_cell[i]))^(n)
        else
            P = Pc
            ρ = ρc
        end
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnP]] = log(P)
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnT]] = log(P * μ / (CGAS * ρ))
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:H1]] = 1.0
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:He4]] = 0
    end

    # set m and dm
    sm.mstar = M
    sm.dm = dms
    sm.m = m_face

    # set luminosity
    for i = 1:(sm.nz - 1)
        μ = 0.5
        Pface = Pc * (θ_n(ξ_face[i]))^(n + 1)
        ρface = ρc * (θ_n(ξ_face[i]))^(n)
        Tface = Pface * μ / (CGAS * ρface)
        dlnT = sm.ind_vars[(i) * sm.nvars + sm.vari[:lnT]] - sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnT]]
        dlnP = sm.ind_vars[(i) * sm.nvars + sm.vari[:lnP]] - sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnP]]
        κ = 0.2
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lum]] = (dlnT / dlnP) *
                                                          (16π * CRAD * CLIGHT * CGRAV * m_face[i] * Tface^4) /
                                                          (3κ * Pface * LSUN)
    end

    # special cases, just copy values at edges
    sm.ind_vars[(sm.nz - 1) * sm.nvars + sm.vari[:lnP]] = sm.ind_vars[(sm.nz - 2) * sm.nvars + sm.vari[:lnP]]
    sm.ind_vars[(sm.nz - 1) * sm.nvars + sm.vari[:lnT]] = sm.ind_vars[(sm.nz - 2) * sm.nvars + sm.vari[:lnT]]
    sm.ind_vars[(sm.nz - 1) * sm.nvars + sm.vari[:lum]] = sm.ind_vars[(sm.nz - 2) * sm.nvars + sm.vari[:lum]]

    sm.time = 0.0
    sm.dt = initial_dt
    sm.model_number = 0
end



"""
    n1_polytrope_initial_condition(sm::StellarModel, M::Real, R::Real; initial_dt=100 * SECYEAR)

Initializes a stellar model `sm` with values corresponding to an n=1 polytrope, setting the independent variables
`sm.ind_vars`, etc. accordingly. Also sets the initial timestep to be taken, `initial_dt`.
"""
function n1_polytrope_initial_condition!(sm::StellarModel, M::Real, R::Real; initial_dt=100 * SECYEAR)
    logdqs = get_logdq.(1:(sm.nz), sm.nz, -3.0, 0.0, 100)
    dqs = 10 .^ logdqs
    dqs = dqs ./ sum(dqs)
    dms = dqs .* M
    m_face = cumsum(dms)
    m_cell = cumsum(dms)
    # correct m_center
    for i = 1:(sm.nz)
        if i == 1
            m_cell[i] = 0
        elseif i != sm.nz
            m_cell[i] = m_cell[i] - 0.5 * dms[i]
        end
    end

    n = 1  # n = 1 polytrope after all...
    rn = R / π  # ξ is defined as r/rn, where rn^2=(n+1)Pc/(4π G ρc^2)

    ρc = M / (4π * rn^3 * (-π^2 * ForwardDiff.derivative(θ_n, π)))
    Pc = 4π * CGRAV * rn^2 * ρc^2 / (n + 1)

    ξ_cell = zeros(sm.nz)
    ξ_face = zeros(sm.nz)
    function mfunc(ξ, m)
        return m - 4π * rn^3 * ρc * (-(-sin(ξ) + cos(ξ) * ξ))
    end

    for i = 1:(sm.nz)
        if i == 1
            ξ_cell[i] = 0
        elseif i == sm.nz
            ξ_cell[i] = π
        else
            mfunc_anon = ξ -> mfunc(ξ, m_cell[i])
            ξ_cell[i] = find_zero(mfunc_anon, (0, π), Bisection())
        end
        if i == sm.nz
            ξ_face[i] = π
        else
            mfunc_anon = ξ -> mfunc(ξ, m_face[i])
            ξ_face[i] = find_zero(mfunc_anon, (0, π), Bisection())
        end
    end

    # set radii, pressure and temperature, assuming ideal gas without Prad
    for i = 1:(sm.nz)
        μ = 0.5
        XH = 1.0
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnr]] = log(rn * ξ_face[i])
        if i > 1
            P = Pc * (θ_n(ξ_cell[i]))^(n + 1)
            ρ = ρc * (θ_n(ξ_cell[i]))^(n)
        else
            P = Pc
            ρ = ρc
        end
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnP]] = log(P)
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnT]] = log(P * μ / (CGAS * ρ))
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:H1]] = 1.0
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:He4]] = 0
    end

    # set m and dm
    sm.mstar = M
    sm.dm = dms
    sm.m = m_face

    # set luminosity
    for i = 1:(sm.nz - 1)
        μ = 0.5
        Pface = Pc * (θ_n(ξ_face[i]))^(n + 1)
        ρface = ρc * (θ_n(ξ_face[i]))^(n)
        Tface = Pface * μ / (CGAS * ρface)
        dlnT = sm.ind_vars[(i) * sm.nvars + sm.vari[:lnT]] - sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnT]]
        dlnP = sm.ind_vars[(i) * sm.nvars + sm.vari[:lnP]] - sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnP]]
        κ = 0.2
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lum]] = (dlnT / dlnP) *
                                                          (16π * CRAD * CLIGHT * CGRAV * m_face[i] * Tface^4) /
                                                          (3κ * Pface * LSUN)
    end

    # special cases, just copy values at edges
    sm.ind_vars[(sm.nz - 1) * sm.nvars + sm.vari[:lnP]] = sm.ind_vars[(sm.nz - 2) * sm.nvars + sm.vari[:lnP]]
    sm.ind_vars[(sm.nz - 1) * sm.nvars + sm.vari[:lnT]] = sm.ind_vars[(sm.nz - 2) * sm.nvars + sm.vari[:lnT]]
    sm.ind_vars[(sm.nz - 1) * sm.nvars + sm.vari[:lum]] = sm.ind_vars[(sm.nz - 2) * sm.nvars + sm.vari[:lum]]

    sm.time = 0.0
    sm.dt = initial_dt
    sm.model_number = 0
end

###########################MY CODE ADDED BELOW##############################################

dydx(x,y,z,n) = z #derivative dy/dx = z
dzdx(x,y,z,n) = -y^n -2*z/x #derivative dz/dx = -y^n -2*z/x
#approximation of y and z for small ξ
y_smallx(x,n) = 1 - 1/6*x^2 + n/120*x^4 -n*(8*n-5)/1520*x^6
z_smallx(x,n) = - 1/3*x + n/30*x^3 -3*n*(8*n-5)/760*x^5;

function endOfLoop!(xvals::LinRange, yvals::Vector{Float64}, zvals::Vector{Float64}, endIndex::Int)
    slope = (xvals[endIndex-1] - yvals[endIndex-2]) / (xvals[endIndex-1] - xvals[endIndex-2])
    xlast = xvals[endIndex-1] - yvals[endIndex-1] / slope
    newxvals = zeros(endIndex)
    newxvals[1:endIndex-1] = xrange[1:endIndex-1]; newxvals[endIndex] = xlast
    return (newxvals, push!(yrange[1:endIndex-1],0.0))
end

function RungeKutta(n::Int64)
    Δx = 1e-4
    nsteps = 200_000 #maximum number of steps
    #initialize first value of y and z using series approximation
    xvals = LinRange(Δx,nsteps*Δx,nsteps)
    yvals = zeros(nsteps); zvals = zeros(nsteps)
    yvals[1] = y_smallx(Δx,n); zvals[1] = z_smallx(Δx,n)
    
    for i in 2:nsteps
        x = xvals[i-1]; y = yvals[i-1]; z = zvals[i-1]
        k₁ = Δx*dydx(x,y,z,n); l₁ = Δx*dzdx(x,y,z,n)
        ynew = y + k₁/2
        if ynew < 0.0
            xvals, yvals = endOfLoop!(xvals,yvals,zvals,i)
            break
        end
        k₂ = Δx*dydx(x+Δx/2,ynew,z+l₁/2,n); l₂ = Δx*dzdx(x+Δx/2,ynew,z+l₁/2,n)
        ynew = y+k₂/2
        if ynew < 0.0
            xvals, yvals = endOfLoop!(xvals,yvals,zvals,i)
            break
        end
        k₃ = Δx*dydx(x+Δx/2,ynew,z+l₂/2,n); l₃ = Δx*dzdx(x+Δx/2,ynew,z+l₂/2,n)
        ynew = y+k₃
        if ynew < 0.0
            xvals, yvals = endOfLoop!(xvals,yvals,zvals,i)
            break
        end
        k₄ = Δx*dydx(x+Δx,ynew,z+l₃,n);l₄ = Δx*dzdx(x+Δx,ynew,z+l₃,n)
        yvals[i] = y+k₁/6+k₂/3+k₃/3+k₄/6 #new y value
        zvals[i] = z+l₁/6+l₂/3+l₃/3+l₄/6 #new z value
    end
    return xvals, yvals
end

function linear_interpolation(xvalues, yvalues)
    function θ_n(x)
        for i in eachindex(xvalues)
            if xvalues[i] <= x < xvalues[i+1]
                return yvalues[i] + (yvalues[i+1] - yvalues[i]) / (xvalues[i+1] - xvalues[i]) * (x - xvalues[i])
            end
        end
    end
    return θ_n
end

function get_θn_and_ξ1(n)
    xvals, yvals = RungeKutta(n)
    return (linear_interpolation(xvalues,yvalues), xvals[end])
end   

function n_polytrope_initial_condition!(sm::StellarModel, M::Real, R::Real; initial_dt=100 * SECYEAR)
    (θ_n, ξ_1) = get_θn_and_ξ1(n) #get θ_n from numerical integration and interpolation
    logdqs = get_logdq.(1:(sm.nz), sm.nz, -3.0, 0.0, 100)
    dqs = 10 .^ logdqs
    dqs = dqs ./ sum(dqs)
    dms = dqs .* M
    m_face = cumsum(dms)
    m_cell = cumsum(dms)
    # correct m_center
    for i = 1:(sm.nz)
        if i == 1
            m_cell[i] = 0
        elseif i != sm.nz
            m_cell[i] = m_cell[i] - 0.5 * dms[i]
        end
    end

    #n = 1  # n = 1 polytrope after all...
    #rn = R / π  # ξ is defined as r/rn, where rn^2=(n+1)Pc/(4π G ρc^2)

    ρc = M / (4π * rn^3 * (-π^2 * ForwardDiff.derivative(θ_n, π)))
    Pc = 4π * CGRAV * rn^2 * ρc^2 / (n + 1)

    ξ_cell = zeros(sm.nz)
    ξ_face = zeros(sm.nz)
    function mfunc(ξ, m)
        return m - 4π * rn^3 * ρc * (-(-sin(ξ) + cos(ξ) * ξ))
    end

    for i = 1:(sm.nz)
        if i == 1
            ξ_cell[i] = 0
        elseif i == sm.nz
            ξ_cell[i] = π
        else
            mfunc_anon = ξ -> mfunc(ξ, m_cell[i])
            ξ_cell[i] = find_zero(mfunc_anon, (0, π), Bisection())
        end
        if i == sm.nz
            ξ_face[i] = π
        else
            mfunc_anon = ξ -> mfunc(ξ, m_face[i])
            ξ_face[i] = find_zero(mfunc_anon, (0, π), Bisection())
        end
    end

    # set radii, pressure and temperature, assuming ideal gas without Prad
    for i = 1:(sm.nz)
        μ = 0.5
        XH = 1.0
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnr]] = log(rn * ξ_face[i])
        if i > 1
            P = Pc * (θ_n(ξ_cell[i]))^(n + 1)
            ρ = ρc * (θ_n(ξ_cell[i]))^(n)
        else
            P = Pc
            ρ = ρc
        end
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnP]] = log(P)
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnT]] = log(P * μ / (CGAS * ρ))
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:H1]] = 1.0
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:He4]] = 0
    end

    # set m and dm
    sm.mstar = M
    sm.dm = dms
    sm.m = m_face

    # set luminosity
    for i = 1:(sm.nz - 1)
        μ = 0.5
        Pface = Pc * (θ_n(ξ_face[i]))^(n + 1)
        ρface = ρc * (θ_n(ξ_face[i]))^(n)
        Tface = Pface * μ / (CGAS * ρface)
        dlnT = sm.ind_vars[(i) * sm.nvars + sm.vari[:lnT]] - sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnT]]
        dlnP = sm.ind_vars[(i) * sm.nvars + sm.vari[:lnP]] - sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnP]]
        κ = 0.2
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lum]] = (dlnT / dlnP) *
                                                          (16π * CRAD * CLIGHT * CGRAV * m_face[i] * Tface^4) /
                                                          (3κ * Pface * LSUN)
    end

    # special cases, just copy values at edges
    sm.ind_vars[(sm.nz - 1) * sm.nvars + sm.vari[:lnP]] = sm.ind_vars[(sm.nz - 2) * sm.nvars + sm.vari[:lnP]]
    sm.ind_vars[(sm.nz - 1) * sm.nvars + sm.vari[:lnT]] = sm.ind_vars[(sm.nz - 2) * sm.nvars + sm.vari[:lnT]]
    sm.ind_vars[(sm.nz - 1) * sm.nvars + sm.vari[:lum]] = sm.ind_vars[(sm.nz - 2) * sm.nvars + sm.vari[:lum]]

    sm.time = 0.0
    sm.dt = initial_dt
    sm.model_number = 0
end

