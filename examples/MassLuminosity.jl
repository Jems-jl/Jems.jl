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
using Jems.StellarModels
using Jems.Evolution
using Jems.ReactionRates


## Plotting with Makie

using CairoMakie, LaTeXStrings, MathTeXEngine
basic_theme = Theme(fonts=(regular=texfont(:text), bold=texfont(:bold),
                           italic=texfont(:italic), bold_italic=texfont(:bolditalic)),
                    fontsize=30, resolution=(1000, 750), linewidth=7,
                    Axis=(xlabelsize=40, ylabelsize=40, titlesize=40, xgridvisible=false, ygridvisible=false,
                          spinewidth=2.5, xminorticksvisible=true, yminorticksvisible=true, xtickalign=1, ytickalign=1,
                          xminortickalign=1, yminortickalign=1, xticksize=14, xtickwidth=2.5, yticksize=14,
                          ytickwidth=2.5, xminorticksize=7, xminortickwidth=2.5, yminorticksize=7, yminortickwidth=2.5,
                          xticklabelsize=35, yticklabelsize=35, xticksmirrored=true, yticksmirrored=true),
                    Legend=(patchsize=(70, 10), framevisible=false, patchlabelgap=20, rowgap=10))
set_theme!(basic_theme)




##

massrange = (0.5:0.5:6)
luminosity = []

for mass in massrange
    println(" ---------------------------------------------------------------------------------")
    println(" ---------------------------------------------------------------------------------")
    println("Mass = $mass")

    ## Model creation

    varnames = [:lnρ, :lnT, :lnr, :lum]
    structure_equations = [Evolution.equationHSE, Evolution.equationT,
                           Evolution.equationContinuity, Evolution.equationLuminosity]
    remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                              StellarModels.split_lnT, StellarModels.split_xa]
    net = NuclearNetwork([:H1,:He4,:C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
    nz = 1000
    nextra = 100
    eos = EOS.IdealEOS(false)
    opacity = Opacity.SimpleElectronScatteringOpacity()
    sm = StellarModel(varnames, structure_equations, nz, nextra,
                      remesh_split_functions, net, eos, opacity);

    ## Initialize StellarModel and evaluate equations and jacobian
    n=3
    StellarModels.n_polytrope_initial_condition!(n, sm, mass*MSUN, 100 * RSUN; initial_dt=10 * SECYEAR)
    Evolution.set_step_info!(sm, sm.esi)
    Evolution.cycle_step_info!(sm);
    Evolution.set_step_info!(sm, sm.ssi)
    Evolution.eval_jacobian_eqs!(sm)

    ## Evolving our model

    open("example_options.toml", "w") do file
        write(file,
              """
              [remesh]
              do_remesh = true

              [solver]
              newton_max_iter_first_step = 1000
              newton_max_iter = 200

              [timestep]
              dt_max_increase = 10.0
              delta_R_limit = 0.01
              delta_Tc_limit = 0.01

              [termination]
              max_model_number = 300
              max_center_T = 4e7

              [plotting]
              do_plotting = false
              wait_at_termination = false
              plotting_interval = 1

              window_specs = ["HR", "profile", "history"]
              window_layouts = [[1, 1],  # arrangement of plots
                                [2, 1],
                                [3, 1]
                                ]

              profile_xaxis = 'mass'
              profile_yaxes = ['log10_T']
              profile_alt_yaxes = ['X','Y']

              history_xaxis = 'star_age'
              history_yaxes = ['R_surf']
              history_alt_yaxes = ['T_center']

              [io]
              profile_interval = 50
              terminal_header_interval = 100
              terminal_info_interval = 100

              """)
    end
    StellarModels.set_options!(sm.opt, "./example_options.toml")
    rm(sm.opt.io.hdf5_history_filename; force=true)
    rm(sm.opt.io.hdf5_profile_filename; force=true)
    StellarModels.n_polytrope_initial_condition!(n, sm, mass*MSUN, 100 * RSUN; initial_dt=1000 * SECYEAR)
    @time sm = Evolution.do_evolution_loop(sm);
    history = StellarModels.get_history_dataframe_from_hdf5("history.hdf5")
    push!(luminosity, history[!, "L_surf"][end])

end
##

f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(M/M_\odot)", ylabel=L"\log_{10}(L/L_\odot)")
scatter!(ax, log10.(massrange), log10.(luminosity))
#linear regression
using LinearAlgebra
A = [ones(length(massrange)) log10.(massrange)]
b = log10.(luminosity)
x = A\b
y = A*x
lines!(ax, log10.(massrange), y, color=:red,linewidth=0.1)
println("Slope = ", x[2])
f
## Plot a funny HR diagram


f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(T_\mathrm{eff}/[K])", ylabel=L"\log_{10}(L/L_\odot)", xreversed=true)
history = StellarModels.get_history_dataframe_from_hdf5("history.hdf5")
lines!(ax, log10.(history[!, "T_surf"]), log10.(history[!, "L_surf"]))
f


##
history[!, "L_surf"][end]

##Cleanup
rm("history.hdf5")
rm("profiles.hdf5")
rm("example_options.toml")
