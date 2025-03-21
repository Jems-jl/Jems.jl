#=
# NuclearBurning.jl

This notebook provides a simple example of a star with simplified microphysics undergoing nuclear burning.
Import all necessary Jems modules. We will also do some benchmarks, so we import BenchmarkTools as well.
=#
using BenchmarkTools
using Jems.Chem
using Jems.Constants
using Jems.EOS
using Jems.Opacity
using Jems.NuclearNetworks
using Jems.Turbulence
using Jems.StellarModels
using Jems.Evolution
using Jems.ReactionRates

##
#=
### Model creation

We start by creating the stellar model. In this example we consider a model with 6 independent variables, two of which
correspond to composition. The independent variables here are $\ln(P)$, $\ln(T)$, $\ln(r)$, the luminosity $L$ and the
mass fractions of Hydrogen and Helium.

The Evolution module has pre-defined equations corresponding to these variables, which we provide here. For now, only a
simple (fully ionized) ideal gas law EOS is available. Similarly, only a simple simple electron scattering opacity equal
to $\kappa=0.2(1+X)\;[\mathrm{cm^2\;g^{-1}}]$ is available.
=#

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
#net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
reactions::Vector{Tuple{Symbol,Symbol}} = []
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16], reactions)
#include("full_net.jl")
nz = 1000
nextra = 100
eos = EOS.IdealEOS(false)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                    nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);

##
#=
### Initialize StellarModel and evaluate equations and jacobian

We do not have a working initial condition yet. We require pressure, temperature profiles. One simple available initial
condition is that of an n=1 polytrope. This sets the pressure and density and computes the temperature from the EOS. The
luminosity is initialized by assuming pure radiative transport for the temperature gradient produced by the polytrope.

The normal evolution loop will store the information at the end of the step into an attribute of type `StellarStepInfo`,
stored at `sm.esi` (_end step info_). After initializing our polytrope we can mimic that behavior by calling 
`set_end_step_info!(sm)`. We then 'cycle' this info into the information of a hypothetical previous step with
`cycle_step_info`, so now `sm.psi` contains our initial condition. Finally we call `set_start_step_info` to use `sm.psi`
(_previous step info_) to populate the information needed before the Newton solver in `sm.ssi` (_start step info_).
At last we are in position to evaluate the equations and compute the Jacobian.
=#
n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
StellarModels.cycle_props!(sm);
StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

##
#=
### Benchmarking

The previous code leaves everything ready to solve the linearized system.
For now we make use of a the serial Thomas algorithm for tridiagonal block matrices.
We first show how long it takes to evaluate the Jacobian matrix. This requires two
steps, the first is to evaluate properties across the model (for example, the EOS)
and then evaluate all differential equations and fill the Jacobian.
=#
@benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
end
##

@benchmark begin
    Evolution.eval_jacobian_eqs!($sm)
end

##
#=

To get an idea of how much a complete iteration of the solver takes, we need to Perform
the jacobian evaluation as a setup for the benchmark. This is because the solver
destroys the Jacobian to perform in-place operations.
=#

@benchmark begin
    Evolution.thomas_algorithm!($sm)
end setup=(Evolution.eval_jacobian_eqs!($sm))

##
#=
### Evolving our model

We can now evolve our star! We will initiate a $1M_\odot$ star with a radius of $100R_\odot$ using an n=1 polytrope (it
would be much better to use n=3 or n=3/2 polytropes, for now I only use this because there is a simple analytical
solution). The star is expected to contract until it ignites hydrogen. We set a few options for the simulation with a
toml file, which we generate dynamically. These simulation should complete in about a thousand steps once it reaches the
`max_center_T` limit.

Output is stored in HDF5 files, and easy to use functions are provided with the Evolution module to turn these HDF5
files into DataFrame objects. HDF5 output is compressed by default.
=#
open("example_options.toml", "w") do file
    write(file,
          """
          [remesh]
          do_remesh = true

          [solver]
          newton_max_iter_first_step = 1000
          initial_model_scale_max_correction = 0.2
          newton_max_iter = 50
          scale_max_correction = 0.1
          report_solver_progress = true
          solver_progress_iter = 1

          [timestep]
          dt_max_increase = 1.5
          dt_retry_decrease = 0.1
          delta_R_limit = 0.01
          delta_Tc_limit = 0.01
          delta_Xc_limit = 0.005

          [termination]
          max_model_number = 3000
          max_center_T = 1e8

          [plotting]
          do_plotting = true
          wait_at_termination = false
          plotting_interval = 1

          window_specs = ["HR", "Kippenhahn", "profile", "TRhoProfile"]
          window_layout = [[1, 1],  # arrangement of plots
                            [1, 2],
                            [2, 1],
                            [2, 2]
                            ]

          profile_xaxis = 'mass'
          profile_yaxes = ['log10_T']
          profile_alt_yaxes = ['X','Y']

          history_xaxis = 'age'
          history_yaxes = ['R_surf']
          history_alt_yaxes = ['T_center']

          max_log_eps = 5.0

          [io]
          profile_interval = 50
          terminal_header_interval = 100
          terminal_info_interval = 100
          hdf5_profile_filename = "1Msun_09.hdf5"

          """)
end
StellarModels.set_options!(sm.opt, "./example_options.toml")
rm(sm.opt.io.hdf5_history_filename; force=true)
rm(sm.opt.io.hdf5_profile_filename; force=true)

n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.9, 0.02, 0.0, Chem.abundance_lists[:ASG_09], 
                                            1 * MSUN, 1 * RSUN; initial_dt=10 * SECYEAR)
@time Evolution.do_evolution_loop!(sm);

##
#=
### Plotting with Makie

Now that our simulation is complete we can analyze the results. We make use of the Makie package for this. I'm not a fan
of the Makie defaults, so I adjust them. I normally also adjust the fonts to be consistent with \LaTeX, but I avoid that
here so we don't need to distribute those fonts together with Jems.
=#
using CairoMakie, LaTeXStrings, MathTeXEngine
basic_theme = Theme(fonts=(regular=texfont(:text), bold=texfont(:bold),
                           italic=texfont(:italic), bold_italic=texfont(:bolditalic)),
                    fontsize=30, size=(1000, 750), linewidth=7,
                    Axis=(xlabelsize=40, ylabelsize=40, titlesize=40, xgridvisible=false, ygridvisible=false,
                          spinewidth=2.5, xminorticksvisible=true, yminorticksvisible=true, xtickalign=1, ytickalign=1,
                          xminortickalign=1, yminortickalign=1, xticksize=14, xtickwidth=2.5, yticksize=14,
                          ytickwidth=2.5, xminorticksize=7, xminortickwidth=2.5, yminorticksize=7, yminortickwidth=2.5,
                          xticklabelsize=35, yticklabelsize=35, xticksmirrored=true, yticksmirrored=true),
                    Legend=(patchsize=(70, 10), framevisible=false, patchlabelgap=20, rowgap=10))
set_theme!(basic_theme)

##
#=
### Compare against polytropes

Below we see how the profile of the star compares to different polytropes. We make use of the facility tools to obtain
DataFrame objects out of the hdf5 output. In particular, `get_profile_names_from_hdf5` will provide the names of all 
profiles contained within the hdf5 file, while `get_profile_dataframe_from_hdf5` is used to obtain one DataFrame
corresponding to one stellar profile. The animation is constructed using the `Observable` type that makie provides. Note
that the zero points of the polytropes are arbitrary.
=#
profile_names = StellarModels.get_profile_names_from_hdf5("profiles.hdf5")

f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(\rho/\mathrm{[g\;cm^{-3}]})", ylabel=L"\log_{10}(P/\mathrm{[dyn]})")

pname = Observable(profile_names[1])

profile = @lift(StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", $pname))
#To see why this is done this way, see https://docs.makie.org/stable/explanations/nodes/index.html#problems_with_synchronous_updates
#the main issue is that remeshing changes the size of the arrays
log10_ρ = @lift($profile[!, "log10_ρ"])
log10_P = Observable(rand(length(log10_ρ.val)))

profile_line = lines!(ax, log10_ρ, log10_P; label="real profile")
xvals = LinRange(-13, 4, 100)
lines!(ax, xvals, (1 + 1 / 1) .* xvals .+ 20; label="n=1")
lines!(ax, xvals, (1 + 1 / (1.5)) .* xvals .+ 15; label="n=1.5")
lines!(ax, xvals, (1 + 1 / 3) .* xvals .+ 15; label="n=3")
axislegend(ax; position=:rb)

model_number_str = @lift("model number=$(parse(Int,$pname))")
profile_text = text!(ax, -10, 20; text=model_number_str)

record(f, "rho_P_evolution.gif", profile_names[1:end]; framerate=2) do profile_name
    profile = StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", profile_name)
    log10_P.val = profile[!, "log10_P"]
    pname[] = profile_name
end

# ![Movie polytrope](./rho_P_evolution.gif)

##
#=
### Check nuclear burning

We see that the structure evolves towards an n=3 polytrope. Deviations near the core are due to the non-homogeneous
composition as hydrogen is burnt. We can similarly visualize how the hydrogen mass fraction changes in the simulation.
In here, only one frame shows the hydrogen that was burnt. To better visualize that you can adjust `profile_interval` in
the [IO](Evolution.md##Io.jl) options (and probably adjust the framerate).
=#
profile_names = StellarModels.get_profile_names_from_hdf5("profiles.hdf5")

f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\mathrm{Mass}\;[M_\odot]", ylabel=L"Abundance")

pname = Observable(profile_names[1])

profile = @lift(StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", $pname))
mass = @lift($profile[!, "mass"])
X = Observable(rand(length(mass.val)))
Y = Observable(rand(length(mass.val)))
model_number_str = @lift("model number=$(parse(Int,$pname))")

profile_line = lines!(ax, mass, X; label="X")
profile_line = lines!(ax, mass, Y; label="Y")
profile_text = text!(ax, 0.7, 0.0; text=model_number_str)
axislegend(ax; position=:rb)

record(f, "X_evolution.gif", profile_names[1:end]; framerate=2) do profile_name
    profile = StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", profile_name)
    X.val = profile[!, "X"]
    Y.val = profile[!, "Y"]
    pname[] = profile_name
end

# ![Movie polytrope](./X_evolution.gif)

##
#=
### Plot a funny HR diagram

Finally, we can also access the history data of the simulation. We use this to plot a simple HR diagram. As our
microphysics are very simplistic, and the initial condition is not very physical, this looks a bit funny!
=#
f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(T_\mathrm{eff}/[K])", ylabel=L"\log_{10}(L/L_\odot)", xreversed=true)
history = StellarModels.get_history_dataframe_from_hdf5("history.hdf5")
lines!(ax, log10.(history[!, "T_surf"]), log10.(history[!, "L_surf"]))
f

##
#=
### Perform some cleanup

Internally we want to prevent storing any of the hdf5 files into our git repos, so I remove them. You can also take
advantage of `julia` as a scripting language to post-process your simulation output in a similar way.
=#
rm("history.hdf5")
rm("profiles.hdf5")
rm("example_options.toml")

##
using Jems.Constants
using Jems.Chem
using Jems.EOS
#function helmholtz_free_energy(lnT::T1, lnρ::T1,
#    xa::AbstractVector{T2}, species::Vector{Symbol},
#    frac_HII::T3, frac_HeII::T3, frac_HeIII::T3) where {T1<:Real, T2<:Real, T3<:Real}
#
#    T = exp(lnT)
#    ρ = exp(lnρ)
#
#    Frad = - CRAD*T^4/(3*ρ) # radiation term 
#    n = 0
#    ne = 0
#    for i in eachindex(species)
#        mj = Chem.isotope_list[species[i]].mass*AMU
#        nj = xa[i]*ρ/mj
#        n += nj
#        if species[i]==:H1
#            ne += frac_HII*nj
#        elseif species[i]==:He4
#            ne += frac_HeII*nj + 2*frac_HeIII*nj
#            #ne += 2*frac_HeIII*nj
#        else
#            ne += nj*Chem.isotope_list[species[i]].Z
#        end
#    end
#    Fideal_ion = 0
#    Fideal_mix = 0
#    mavg = 0
#    for i in eachindex(species)
#        mj = Chem.isotope_list[species[i]].mass*AMU
#        nj = xa[i]*ρ/mj
#        yj = nj/n
#        mavg += yj*mj
#
#        nQj = (2*π*Constants.HBAR^2/(mj*KERG*T))^(-3/2)
#
#        Fideal_ion += yj*(log((nj+1e-99)/nQj)-1)
#
#        Fideal_mix += yj*log(yj+1e99)
#    end
#    Fideal_ion = KERG*T/mavg*Fideal_ion
#    Fideal_mix = KERG*T/mavg*Fideal_mix
#
#    nQj = (2*π*Constants.HBAR^2/(ME*KERG*T))^(-3/2)
#    Felectron = ne*KERG*T/ρ*(log((ne+1e-99)/nQj)-1)
#
#    frac_HI = 1 - frac_HII
#    frac_HeI = 1 - frac_HeII - frac_HeIII 
#    # do HI
#    mH = Chem.isotope_list[:H1].mass*AMU
#    iH1 = findfirst(==(:H1), species)
#    nHI = frac_HI*xa[iH1]*ρ/mH
#    Fideal_ion += -nHI*13.6*Constants.EV2ERG/ρ
#    # do HeI
#    mHe = Chem.isotope_list[:He4].mass*AMU
#    iHe4 = findfirst(==(:He4), species)
#    nHeI = frac_HeI*xa[iHe4]*ρ/mHe
#    Fideal_ion += -nHeI*78.98*Constants.EV2ERG/ρ
#    ## do HeII
#    mHe = Chem.isotope_list[:He4].mass*AMU
#    iHe4 = findfirst(==(:He4), species)
#    nHeII = frac_HeII*xa[iHe4]*ρ/mHe
#    Fideal_ion += -nHeII*54.40*Constants.EV2ERG/ρ
#
#    return Frad + Fideal_ion + Fideal_mix + Felectron
#end
function helmholtz_free_energy(lnT::T1, lnρ::T1,
    xa::AbstractVector{T2}, species::Vector{Symbol},
    frac_HII::T3, frac_HeII::T3, frac_HeIII::T3) where {T1<:Real, T2<:Real, T3<:Real}

    T = exp(lnT)
    ρ = exp(lnρ)

    frac_HI = max(0.0,min(1.0,1 - frac_HII))
    frac_HeI = max(0.0, min(1.0,1 - frac_HeII - frac_HeIII)) 

    Frad = - CRAD*T^4/(3*ρ) # radiation term 

    ne = 0
    Fideal_ion = 0
    for i in eachindex(species)
        mj = Chem.isotope_list[species[i]].mass*AMU
        nj = xa[i]/mj

        if species[i]==:H1
            Fideal_ion += frac_HI*nj*(log((ρ*(frac_HI*nj+1e-99)*Constants.PLANCK_H^3)/(2*2*pi*mj*KERG*T)^(3/2))-1-13.6*Constants.EV2ERG/(KERG*T))
            Fideal_ion += frac_HII*nj*(log((ρ*(frac_HII*nj+1e-99)*Constants.PLANCK_H^3)/(2*pi*mj*KERG*T)^(3/2))-1)
            ne += frac_HII*nj
        elseif species[i]==:He4
            Fideal_ion += frac_HeI*nj*(log((ρ*(frac_HeI*nj+1e-99)*Constants.PLANCK_H^3)/(2*pi*mj*KERG*T)^(3/2))-1-78.98*Constants.EV2ERG/(KERG*T))
            Fideal_ion += frac_HeII*nj*(log((ρ*(frac_HeII*nj+1e-99)*Constants.PLANCK_H^3)/(2*2*pi*mj*KERG*T)^(3/2))-1-54.40*Constants.EV2ERG/(KERG*T))
            Fideal_ion += frac_HeIII*nj*(log((ρ*(frac_HeIII*nj+1e-99)*Constants.PLANCK_H^3)/(2*pi*mj*KERG*T)^(3/2))-1)
            ne += frac_HeII*nj + 2*frac_HeIII*nj
        else
            Fideal_ion += nj*(log((ρ*(nj+1e-99)*PLANCK_H^3)/(2*pi*mj*KERG*T)^(3/2))-1)
            ne += nj*Chem.isotope_list[species[i]].Z
        end
    end
    Fideal_ion = KERG*T*Fideal_ion

    Felectron = KERG*T*ne*(log((ρ*(ne+1e-99)*Constants.PLANCK_H^3)/(2*pi*ME*KERG*T)^(3/2))-1)

    return Frad + Fideal_ion + Felectron
end

using Optim, LineSearches
function ion_fracs(lnT::TT, lnρ::TT,
    xa::AbstractVector{TTT}, species::Vector{Symbol}) where{TT, TTT}
    f(x) = helmholtz_free_energy(lnT, lnρ, xa, species, x[1], (1-x[3])*x[2], x[3])
    frac_HII_guess::TT = 0.5
    frac_HeII_guess::TT = 0.5
    frac_HeIII_guess::TT = 0.5

    lower = [0.00001,0.00001,0.00001]
    upper = [0.99999,0.99999,0.99999]
    inner_optimizer = GradientDescent()
    res = optimize(f, lower, upper, [frac_HII_guess, frac_HeII_guess, frac_HeIII_guess],
                    Fminbox(inner_optimizer))
    vals = Optim.minimizer(res)
    return [vals[1], (1-vals[3])*vals[2], vals[3]]
end

function minimized_helmholtz_free_energy(lnT::TT, lnρ::TT,
    xa::AbstractVector{TTT}, species::Vector{Symbol}) where{TT, TTT}
    fracs = ion_fracs(lnT, lnρ, xa, species)
    return helmholtz_free_energy(lnT, lnρ, xa, species, fracs[1], fracs[2], fracs[3])
end
##
using ForwardDiff
xa = [0.7, 0.3]
species = [:H1, :He4]
T = 1e6
ρ = 1e-6
f = x->minimized_helmholtz_free_energy(log(x[1]), log(x[2]), xa, species)
grad = ForwardDiff.gradient(f, [T,ρ])
hessian = ForwardDiff.hessian(f, [T,ρ])

##
F = minimized_helmholtz_free_energy(log(T), log(ρ), xa, species)
P = ρ^2*grad[2]
s = -grad[1]
e = F + T*s
cv = grad[1]+s-T*hessian[1,1]
χT = T*ρ^2/P*hessian[2,1]
χρ = ρ/P*(2*ρ*grad[2]+ρ^2*hessian[2,2])
Γ3 = 1 + P/(ρ*cv*T)*χT
Γ1 = χρ + (Γ3 - 1)*χT
∇ad = (Γ3 - 1)/Γ1

##
nvals = 50
logT_vals = LinRange(3,4,nvals)
∇_Helm = zeros(nvals)
∇_ion = zeros(nvals)

eos = EOS.IdealEOS(true)
eos_res = EOSResults{Float64}()

ρ = 1e-8
xa = [0.7, 0.3]
species = [:H1, :He4]
for (i,logT) in enumerate(logT_vals)
    T = 10^(logT)
    f = x->minimized_helmholtz_free_energy(log(x[1]), log(x[2]), xa, species)
    grad = ForwardDiff.gradient(f, [T,ρ])
    hessian = ForwardDiff.hessian(f, [T,ρ])
    F = minimized_helmholtz_free_energy(log(T), log(ρ), xa, species)
    P = ρ^2*grad[2]
    s = -grad[1]
    e = F + T*s
    cv = grad[1]+s-T*hessian[1,1]
    χT = T*ρ^2/P*hessian[2,1]
    χρ = ρ/P*(2*ρ*grad[2]+ρ^2*hessian[2,2])
    Γ3 = 1 + P/(ρ*cv*T)*χT
    Γ1 = χρ + (Γ3 - 1)*χT
    ∇ad = (Γ3 - 1)/Γ1

    set_EOS_resultsTρ!(eos, eos_res, log(T), log(ρ), xa, species)

    ∇_Helm[i] = ∇ad
    ∇_ion[i] = eos_res.∇ₐ
end
##
using CairoMakie
f = Figure()
ax = Axis(f[1,1])
lines!(logT_vals, ∇_ion)
lines!(logT_vals, ∇_Helm)

f

##
eos = EOS.IdealEOS(true)
eos_res = EOSResults{Float64}()
set_EOS_resultsTρ!(eos, eos_res, log(T), log(ρ), xa, species)
eos_res.P

#using ForwardDiff
#using LinearAlgebra
#function ion_fracs(lnT::TT, lnρ::TT,
#    xa::AbstractVector{TTT}, species::Vector{Symbol}) where{TT, TTT}
#    f(x) = helmholtz_free_energy(lnT, lnρ, xa, species, x[1], (1-x[3])*x[2], x[3])
#    guess = ones(TT,3)*0.01
#    guess[2] = 0.0
#
#    new_guess = ones(TT,3)*0.5
#
#    for i in 1:100
#        F = f(guess)
#        grad = ForwardDiff.gradient(f,guess)
#        hessian = ForwardDiff.hessian(f,guess)
#        corr = -hessian\grad
#        new_guess = guess + corr
#        for i in 1:3
#            if new_guess[i] < 0
#                new_guess[i] = 0
#            elseif new_guess[i] > 1
#                new_guess[i] = 1
#            end
#        end
#        eff_corr = new_guess - guess
#        corr_mag = norm(eff_corr)
#        if corr_mag > 0.05
#            eff_corr = eff_corr*0.05/corr_mag
#        end
#        guess = guess + eff_corr
#
#        #@show guess, F
#    end
#    return [guess[1], (1-guess[3])*guess[2], guess[3]]
#end

##
using ForwardDiff
xa = [0.7, 0.3];species=[:H1,:He4]
Tnormal = 1e5
ρnormal = 1e-10
T = ForwardDiff.Dual{:caca}(Tnormal,ForwardDiff.Partials((1.0, 0.0)))
ρ = ForwardDiff.Dual{:caca}(ρnormal ,ForwardDiff.Partials((0.0, 1.0)))
F =ion_fracs(log(Tnormal), log(ρnormal), xa, species)
##
helmholtz_free_energy(log(Tnormal), log(ρnormal), xa, species, 1.0, 0.0, 1.0)

##
using CairoMakie
f = Figure()
ax = Axis(f[1,1], xlabel="logrho", ylabel="logT")
logρ_vals = LinRange(-10,5,50)
logT_vals = LinRange(3,6,50)
for logρ in logρ_vals
    for logT in logT_vals
        fracs = ion_fracs(logT*log(10), logρ*log(10), xa, species)
        if fracs[1] > 0.05 && fracs[1] < 0.95
            scatter!(ax, [logρ], [logT], color="red")
        end
        if fracs[2] > 0.05 && fracs[2] < 0.95
            scatter!(ax, [logρ], [logT], color="blue")
        end
        if fracs[3] > 0.05 && fracs[3] < 0.95
            scatter!(ax, [logρ], [logT], color="gray")
        end
    end
end

save("ionization_regions.png",f)

f
##
using CairoMakie
f = Figure()
ax = Axis(f[1,1], xlabel="logT", ylabel="fraction", title="X=0.7, Y=0.3, rho=1e-10")
logρ = -10
nvals = 1000
logT_vals = LinRange(3,5,nvals)
fracs_HII = zeros(nvals)
fracs_HeII = zeros(nvals)
fracs_HeIII = zeros(nvals)

for i in eachindex(logT_vals)
    logT = logT_vals[i]

    fracs = ion_fracs(logT*log(10), logρ*log(10), xa, species)
    fracs_HII[i] = fracs[1]
    fracs_HeII[i] = fracs[2]
    fracs_HeIII[i] = fracs[3]

end

lines!(ax, logT_vals, fracs_HII, label="HII")
lines!(ax, logT_vals, fracs_HeII, label="HeII") 
lines!(ax, logT_vals, fracs_HeIII, label="HeIII")

save("ionization_fractions.png",f)

f

##
using ForwardDiff
xa = [1.0, 0.0];species=[:H1,:He4]
Tnormal = 1e3
ρnormal = 5e-10
T = ForwardDiff.Dual(Tnormal,ForwardDiff.Partials((1.0, 0.0)))
ρ = ForwardDiff.Dual(ρnormal ,ForwardDiff.Partials((0.0, 1.0)))
F =helmholtz_free_energy(log(T), log(ρ), xa, species, 1.0, 0.0, 1.0)

(ρ^2*F.partials[2]).value

##
eos = EOS.IdealEOS(true)
eos_res = EOSResults{Float64}()
set_EOS_resultsTρ!(eos, eos_res, log(Tnormal), log(ρnormal ), xa, species)
eos_res.P


##
x = ForwardDiff.Dual(2,
        ForwardDiff.Partials((ForwardDiff.Dual(1.0, ForwardDiff.Partials((0.0,0.0))),
                              ForwardDiff.Dual(0.0, ForwardDiff.Partials((0.0,0.0)))))
    )

x^3