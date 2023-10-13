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
    logdqs = zeros(sm.nz + sm.nextra)
    for i in 1:sm.nz
        logdqs[i] = get_logdq(i, sm.nz, -5.0, 0.0, 100)
    end
    dqs = 10 .^ logdqs
    dqs[sm.nz+1:end] .= 0 # extra entries beyond nz have no mass
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
    mfunc_anon = ξ -> mfunc(ξ, 0.99999*M)

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
