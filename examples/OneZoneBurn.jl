#=
# OneZoneBurn.jl

This notebook provides a simple example of a single zone undergoing nuclear burning.
Just as in NuclearBurning.jl, we start by importing all necessary Jems modules.
=#
using Jems.NuclearNetworks
using Jems.StellarModels
using Jems.Evolution
using Jems.Constants
using Jems.DualSupport
using Jems.Chem
using Jems.Plotting
using BenchmarkTools

##
#=
### Model creation

We start by creating our OneZone model. Contrary to a fully-fledged StellarModel, we need only to specify the nuclear
net used, along with its reactions. We also provide a custom equation for the composition, that does not include mixing
(we don't need it, of course, there is only one zone in this model).
=#
net = merge_nuclear_networks([NuclearNetworks.networks[:JINA_PPI],
                              NuclearNetworks.networks[:JINA_PPII],
                              NuclearNetworks.networks[:JINA_PPIII],
                              NuclearNetworks.networks[:JINA_PPIV],
                              NuclearNetworks.networks[:JINA_CNOI],
                              NuclearNetworks.networks[:JINA_CNOII],
                              ])

oz = OneZone(DefaultOneZoneEquationSet(), net);

##
#=
### Set initial conditions
We now provide the initial conditions of our One Zone model, setting the temperature, density, starting dt and age, and
the initial abundances for all isotopes defined in our network above.
=#
oz.props.T = 4e7  # K
oz.props.ρ = 1  # g cm^{-3}
oz.props.ind_vars = zeros(oz.network.nspecies)
# oz.props.ind_vars[oz.network.xa_index[:H1]] = 1.0 # for full hydrogen mixture
mass_fractions = get_mass_fractions(Chem.abundance_lists[:ASG_09],
                                        oz.network.species_names, 0.7, 0.02, 0.0) 
for name in oz.network.species_names
    oz.props.ind_vars[oz.network.xa_index[name]] = mass_fractions[name]
end

oz.props.dt_next = 1 * SECYEAR
oz.props.time = 0.0
oz.props.model_number = 0

open("example_options.toml", "w") do file
    write(file,
          """

          [solver]
          newton_max_iter_first_step = 1000
          initial_model_scale_max_correction = 0.2
          newton_max_iter = 50
          scale_max_correction = 0.2
          report_solver_progress = false
          solver_progress_iter = 50

          [timestep]
          dt_max_increase = 1.5
          delta_Xc_limit = 0.005

          [termination]
          max_model_number = 200

          [io]
          profile_interval = 50
          terminal_header_interval = 100
          terminal_info_interval = 100
          history_values = ["age", "dt", "model_number", "T", "rho",
                            "H1", "D2", "He3", "He4", "Li7", "Be7", "C12", "N14", "O16"]

          """)
end
StellarModels.set_options!(oz.opt, "./example_options.toml")
rm(oz.opt.io.hdf5_history_filename; force=true)

using GLMakie
set_theme!(Plotting.basic_theme())
f = Figure(size=(1400,750))
plots = [Plotting.HistoryPlot(f[1,1], oz, x_name="age", y_name="H1", othery_name="He4", link_yaxes=true)]
plotter = Plotting.Plotter(fig=f,plots=plots)

Evolution.do_one_zone_burn!(oz, plotter=plotter)

##
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
# GLMakie.activate!()
##
### Plot the history
f = Figure();
ax = Axis(f[1, 1]; xlabel="age (year)", ylabel=L"\log_{10}(X)", xscale=log10, xminorticks=IntervalsBetween(10, false))
history = StellarModels.get_history_dataframe_from_hdf5("history.hdf5")
lines!(ax, history[!, "age"], log10.(max.(1e-99,history[!, "H1"])), label=L"^1H")
lines!(ax, history[!, "age"], log10.(max.(1e-99,history[!, "D2"])), label=L"^2H")
lines!(ax, history[!, "age"], log10.(max.(1e-99,history[!, "He3"])), label=L"^3He")
lines!(ax, history[!, "age"], log10.(max.(1e-99,history[!, "He4"])), label=L"^4He")
lines!(ax, history[!, "age"], log10.(max.(1e-99,history[!, "Li7"])), label=L"^7Li")
lines!(ax, history[!, "age"], log10.(max.(1e-99,history[!, "Be7"])), label=L"^7Be")
lines!(ax, history[!, "age"], log10.(max.(1e-99,history[!, "C12"])), label=L"^{12}C")
lines!(ax, history[!, "age"], log10.(max.(1e-99,history[!, "N14"])), label=L"^{14}N")
lines!(ax, history[!, "age"], log10.(max.(1e-99,history[!, "O16"])), label=L"^{16}O")
axislegend(position=:lt)
ylims!(ax, -5,0.1)
xlims!(ax, 1, 1e7)
save("abundance_evolution.png", f)
f

##
#=
### Perform some cleanup

Internally we want to prevent storing any of the hdf5 files into our git repos, so I remove them. You can also take
advantage of `julia` as a scripting language to post-process your simulation output in a similar way.
=#
rm("history.hdf5")
rm("example_options.toml")