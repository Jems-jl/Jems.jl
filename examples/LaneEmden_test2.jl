using BenchmarkTools
using Jems.Chem
using Jems.Constants
using Jems.EOS
using Jems.Opacity
using Jems.NuclearNetworks
using Jems.StellarModels
using Jems.Evolution
using Jems.ReactionRates
using Interpolations
###My code


function RungeKutta(n)

    dydx(x,y,z,n) = z 
    dzdx(x,y,z,n) = -y^n -2*z/x 
    y_smallx(x,n) = 1 - 1/6*x^2 + n/120*x^4 -n*(8*n-5)/1520*x^6
    z_smallx(x,n) = - 1/3*x + n/30*x^3 -3*n*(8*n-5)/760*x^5;

    function endOfLoop!(xvals::LinRange, yvals::Vector{Float64}, zvals::Vector{Float64}, endIndex::Int)
        slope = (yvals[endIndex-1] - yvals[endIndex-2]) / (xvals[endIndex-1] - xvals[endIndex-2])
        xlast = xvals[endIndex-1] - yvals[endIndex-1] / slope
        newxvals = zeros(endIndex)
        newxvals[1:endIndex-1] = xvals[1:endIndex-1]; newxvals[endIndex] = xlast
        #add last entry
        yvals[endIndex] = 0.0
        zvals[endIndex] = zvals[endIndex-1]
        #put first entry (core boundary conditions)
        pushfirst!(yvals,1.0)
        pushfirst!(newxvals,0.0)
        pushfirst!(zvals,0.0)
        return (newxvals,yvals[1:endIndex+1],zvals[1:endIndex+1])
    end

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
            xvals, yvals, zvals = endOfLoop!(xvals,yvals,zvals,i)
            break
        end
        k₂ = Δx*dydx(x+Δx/2,ynew,z+l₁/2,n); l₂ = Δx*dzdx(x+Δx/2,ynew,z+l₁/2,n)
        ynew = y+k₂/2
        if ynew < 0.0
            xvals, yvals, zvals = endOfLoop!(xvals,yvals,zvals,i)
            break
        end
        k₃ = Δx*dydx(x+Δx/2,ynew,z+l₂/2,n); l₃ = Δx*dzdx(x+Δx/2,ynew,z+l₂/2,n)
        ynew = y+k₃
        if ynew < 0.0
            xvals, yvals, zvals = endOfLoop!(xvals,yvals,zvals,i)
            break
        end
        k₄ = Δx*dydx(x+Δx,ynew,z+l₃,n);l₄ = Δx*dzdx(x+Δx,ynew,z+l₃,n)
        yvals[i] = y+k₁/6+k₂/3+k₃/3+k₄/6 #new y value
        zvals[i] = z+l₁/6+l₂/3+l₃/3+l₄/6 #new z value
        #trycatch #still to do
    end
    return xvals, yvals, zvals
end

#function linear_interpolation(xvalues, yvalues)
#    function θ_n(x)
#        for i in 1:length(xvalues)-1
#            if xvalues[i] <= x <= xvalues[i+1]
#                return yvalues[i] + (yvalues[i+1] - yvalues[i]) / (xvalues[i+1] - xvalues[i]) * (x - xvalues[i])
#            end
#        end
#    end
#    return θ_n
#end



function get_logdq(k::Int, nz::Int, logdq_center::TT, logdq_mid::TT, logdq_surf::TT, numregion::Int)::TT where {TT<:Real}
    if k <= numregion
        return logdq_center + (k - 1) * (logdq_mid - logdq_center) / numregion
    elseif k < nz - numregion
        return logdq_mid
    else
        return logdq_mid + (logdq_surf - logdq_mid) * (k - (nz - numregion)) / numregion
    end
end

function n_polytrope_initial_condition!(sm::StellarModel, M::Real, R::Real; initial_dt=100 * SECYEAR)
    n=1 #remove this 
    xvals, yvals, zvals = RungeKutta(n)
    (θ_n, ξ_1, derivative_θ_n) = (linear_interpolation(xvals,yvals), xvals[end],linear_interpolation(xvals,zvals))
    logdqs = zeros(length(sm.dm))
    for i in 1:sm.nz
        logdqs[i] = get_logdq(i, sm.nz, -10.0, 0.0, -6.0, 200)
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

    
    rn = R / ξ_1 # ξ is defined as r/rn, where rn^2=(n+1)Pc/(4π G ρc^2)

    #ρc = M / (4π * rn^3 * (-π^2 * ForwardDiff.derivative(θ_n, π)))
    ρc = M / (4π * rn^3 * (-ξ_1^2 * derivative_θ_n(ξ_1)))
    Pc = 4π * CGRAV * rn^2 * ρc^2 / (n + 1)

    ξ_cell = zeros(sm.nz)
    ξ_face = zeros(sm.nz)
    function mfunc(ξ, m)
        return m - 4π * rn^3 * ρc * (-ξ^2 * derivative_θ_n(ξ))
    end

    for i = 1:(sm.nz)
        if i == 1
            ξ_cell[i] = 0
        elseif i == sm.nz
            ξ_cell[i] = ξ_1
        else
            mfunc_anon = ξ -> mfunc(ξ, m_cell[i])
            ξ_cell[i] = find_zero(mfunc_anon, (0, ξ_1), Bisection())
        end
        if i == sm.nz
            ξ_face[i] = ξ_1
        else
            mfunc_anon = ξ -> mfunc(ξ, m_face[i])
            ξ_face[i] = find_zero(mfunc_anon, (0, ξ_1), Bisection())
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
        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnρ]] = log(ρ)
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
        if i != 1
            dlnP = log(Pc * (θ_n(ξ_cell[i+1]))^(n + 1)) - log(Pc * (θ_n(ξ_cell[i]))^(n + 1))
        else
            dlnP = log(Pc * (θ_n(ξ_cell[i+1]))^(n + 1)) - log(Pc)
        end
        κ = 0.2

        sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lum]] = (dlnT / dlnP) *
                                                          (16π * CRAD * CLIGHT * CGRAV * m_face[i] * Tface^4) /
                                                          (3κ * Pface * LSUN)
    end

    # special cases, just copy values at edges
    sm.ind_vars[(sm.nz - 1) * sm.nvars + sm.vari[:lnρ]] = sm.ind_vars[(sm.nz - 2) * sm.nvars + sm.vari[:lnρ]]
    sm.ind_vars[(sm.nz - 1) * sm.nvars + sm.vari[:lnT]] = sm.ind_vars[(sm.nz - 2) * sm.nvars + sm.vari[:lnT]]
    #sm.ind_vars[(sm.nz - 1) * sm.nvars + sm.vari[:lum]] = sm.ind_vars[(sm.nz - 2) * sm.nvars + sm.vari[:lum]]
    sm.ind_vars[(sm.nz - 1) * sm.nvars + sm.vari[:lum]] = sm.ind_vars[(sm.nz - 3) * sm.nvars + sm.vari[:lum]]
    sm.ind_vars[(sm.nz - 2) * sm.nvars + sm.vari[:lum]] = sm.ind_vars[(sm.nz - 3) * sm.nvars + sm.vari[:lum]]

    sm.time = 0.0
    sm.dt = initial_dt
    sm.model_number = 0
end

###Testing My code using some stuff from NuclearBurngin.jl
using BenchmarkTools
using Jems.Chem
using Jems.Constants
using Jems.EOS
using Jems.Opacity
using Jems.NuclearNetworks
using Jems.StellarModels
using Jems.Evolution
using Jems.ReactionRates
##
varnames = [:lnρ, :lnT, :lnr, :lum]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1,:He4], [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(false)
opacity = Opacity.SimpleElectronScatteringOpacity()
sm = StellarModel(varnames, structure_equations, nz, nextra,
                  remesh_split_functions, net, eos, opacity);

n_polytrope_initial_condition!(sm, MSUN, 100 * RSUN; initial_dt=10 * SECYEAR)
#StellarModels.n1_polytrope_initial_condition!(sm, MSUN, 100 * RSUN; initial_dt=10 * SECYEAR)
Evolution.set_step_info!(sm, sm.esi)
Evolution.cycle_step_info!(sm);
Evolution.set_step_info!(sm, sm.ssi)
Evolution.eval_jacobian_eqs!(sm)
##
sm.ssi.lnP[1:1000]
sm.ssi.lnρ[1:1000]
sm.ssi.L[1:1000]
sm.ssi.lnT[1:1000]
sm.dm[1:1000]

#[sm.ind_vars[sm.nvars*(i-1)+4] for i in 1:nz]


##
using ForwardDiff
using Roots

##
n=3
xvals, yvals, zvals = RungeKutta(n)
(θ_n, ξ_1, derivative_θ_n) = (linear_interpolation(xvals,yvals), xvals[end],linear_interpolation(xvals,zvals))
##
xvals
##
open("example_options.toml", "w") do file
    write(file,
          """
          [remesh]
          do_remesh = true

          [solver]
          newton_max_iter_first_step = 1000
          newton_max_iter = 200

          [timestep]
          dt_max_increase = 2.0

          [termination]
          max_model_number = 2000
          max_center_T = 4e7

          [io]
          profile_interval = 50
          """)
end
StellarModels.set_options!(sm.opt, "./example_options.toml")
rm(sm.opt.io.hdf5_history_filename; force=true)
rm(sm.opt.io.hdf5_profile_filename; force=true)
n_polytrope_initial_condition!(sm, 1*MSUN, 100 * RSUN; initial_dt=1000 * SECYEAR)
@time sm = Evolution.do_evolution_loop(sm);







using CairoMakie
f=Figure()
ax = Axis(f[1,1])
#lines!(ax,xvals, zvals)
lines!(ax, xvals, derivative_θ_n.(xvals))
lines!(ax,xvals, yvals)
lines!(ax,xvals, θ_n.(xvals))

f
##
println(θ_n(0.00001))
##
xvals[1]