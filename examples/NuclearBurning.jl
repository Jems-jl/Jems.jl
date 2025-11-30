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
using Jems.Plotting

##
#=
### Model creation

We start by creating the stellar model. In this example we consider a model with 6 independent variables, two of which
correspond to composition. The independent variables here are $\ln(P)$, $\ln(T)$, $\ln(r)$, the luminosity $L$ and the
mass fractions of Hydrogen and Helium.

The Evolution module has pre-defined equations corresponding to these variables, which we provide here. For now, only a
simple (fully ionized) ideal gas law EOS is available. Similarly, only a simple simple electron scattering opacity equal
to $\kappa=0.2(1+X)\;[\text{cm^2\;g^{-1}}]$ is available.
=#

##

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
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
Information of the model at its present and following step are required at the beginning, the function
`compute_starting_model_properties!` takes care of setting this up.
=#
n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
Evolution.compute_starting_model_properties!(sm)

##
#=
### Benchmarking

The previous code leaves everything ready to solve the linearized system.
For now we make use of a the serial Thomas algorithm for tridiagonal block matrices.
We first show how long it takes to evaluate the Jacobian matrix. This requires two
steps, the first is to evaluate properties across the model (for example, the EOS)
and then evaluate all differential equations and fill the Jacobian. We first benchmark
the evaluation of model properties:
=#
@benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
end

##
#=
And next we benchmark the evaluation of the model equations and construction of the Jacobian:
=#
@benchmark begin
    Evolution.eval_jacobian_eqs!($sm)
end

##
#=

To benchmark the linear solver itself we need to perform
the jacobian evaluation as a setup for the benchmark. This is because the solver
destroys the Jacobian to perform in-place operations.
=#

@benchmark begin
    Evolution.thomas_algorithm!($sm)
end setup=(Evolution.eval_jacobian_eqs!($sm))

##
#=
### Evolving our model

We can now evolve our star! We will initiate a $1M_\odot$ star with a radius of $100R_\odot$ using an n=3 polytrope.
The star is expected to contract until it ignites hydrogen. We set a few options for the simulation with a
toml file, which we generate dynamically. These simulation should complete in about a thousand steps once it reaches the
`max_center_T` limit.

Output is stored in HDF5 files, and easy to use functions are provided with the `StellarModels` module to turn these HDF5
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
          newton_max_iter = 10
          scale_max_correction = 0.1

          [timestep]
          dt_max_increase = 1.5
          delta_R_limit = 0.01
          delta_Tc_limit = 0.01
          delta_Xc_limit = 0.005

          [termination]
          max_model_number = 2000
          max_center_T = 1e8

          [io]
          profile_interval = 50
          terminal_header_interval = 100
          terminal_info_interval = 100

          """)
end
StellarModels.set_options!(sm.opt, "./example_options.toml")
rm(sm.opt.io.hdf5_history_filename; force=true)
rm(sm.opt.io.hdf5_profile_filename; force=true)

#Configure live plots. To turn off one can use `plotter = Plotting.NullPlotter()`
using GLMakie
set_theme!(Plotting.basic_theme())
f = Figure(size=(1400,750))
plots = [Plotting.HRPlot(f[1,1]),
         Plotting.TRhoProfile(f[1,2]),
         Plotting.KippenLine(f[2,1], xaxis=:time, time_units=:Gyr),
         Plotting.AbundancePlot(f[2,2],net,log_yscale=true, ymin=1e-3),
         Plotting.HistoryPlot(f[1,3], sm, x_name="age", y_name="X_center", othery_name="Y_center", link_yaxes=true),
         Plotting.ProfilePlot(f[2,3], sm, x_name="mass", y_name="log10_rho", othery_name="log10_T")]
plotter = Plotting.Plotter(fig=f,plots=plots)

#set initial condition and run model
n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], 
                                            1 * MSUN, 100 * RSUN; initial_dt=10 * SECYEAR)
@time Evolution.do_evolution_loop!(sm, plotter=plotter);

##
#=
### Plotting with Makie

Now that our simulation is complete we can analyze the results. We make use of the Makie package for this. I'm not a fan
of the Makie defaults, so I adjust them.
=#
using CairoMakie, LaTeXStrings
set_theme!(Plotting.basic_theme())

##
#=
#### Compare against polytropes

Below we see how the profile of the star compares to different polytropes. We make use of the facility tools to obtain
DataFrame objects out of the hdf5 output. In particular, `get_profile_names_from_hdf5` will provide the names of all 
profiles contained within the hdf5 file, while `get_profile_dataframe_from_hdf5` is used to obtain one DataFrame
corresponding to one stellar profile. The animation is constructed using the `Observable` type that makie provides. Note
that the zero points of the polytropes are arbitrary.
=#
profile_names = StellarModels.get_profile_names_from_hdf5("profiles.hdf5")

f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(\rho/\text{[g\;cm^{-3}]})", ylabel=L"\log_{10}(P/\text{[dyn]})")

pname = Observable(profile_names[1])

profile = @lift(StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", $pname))
log10_ρ = @lift($profile[!, "log10_rho"])
log10_P = @lift($profile[!, "log10_P"])

profile_line = lines!(ax, log10_ρ, log10_P; label="real profile")
xvals = LinRange(-13, 4, 100)
lines!(ax, xvals, (1 + 1 / 1) .* xvals .+ 20; label="n=1")
lines!(ax, xvals, (1 + 1 / (1.5)) .* xvals .+ 15; label="n=1.5")
lines!(ax, xvals, (1 + 1 / 3) .* xvals .+ 15; label="n=3")
axislegend(ax; position=:rb)

model_number_str = @lift("model number=$(parse(Int,$pname))")
profile_text = text!(ax, -10, 20; text=model_number_str)

record(f, "rho_P_evolution.gif", profile_names[1:end]; framerate=4) do profile_name
    pname[] = profile_name
end

# ![Movie polytrope](./rho_P_evolution.gif)

##
#=
#### Check nuclear burning

We see that the structure evolves towards an n=3 polytrope. Deviations near the core are due to the non-homogeneous
composition as hydrogen is burnt. We can similarly visualize how the hydrogen mass fraction changes in the simulation.
In here, only one frame shows the hydrogen that was burnt. To better visualize that you can adjust `profile_interval` in
the [IO](Evolution.md##Io.jl) options (and probably adjust the framerate).
=#
profile_names = StellarModels.get_profile_names_from_hdf5("profiles.hdf5")

f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\text{Mass}\;[M_\odot]", ylabel=L"\text{Abundance}")

pname = Observable(profile_names[1])

profile = @lift(StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", $pname))
mass = @lift($profile[!, "mass"])
X = @lift($profile[!, "X"])
Y = @lift($profile[!, "Y"])
model_number_str = @lift("model number=$(parse(Int,$pname))")

profile_line = lines!(ax, mass, X; label="X")
profile_line = lines!(ax, mass, Y; label="Y")
profile_text = text!(ax, 0.7, 0.95; text=model_number_str)
axislegend(ax; position=:rb)
ylims!(ax,-0.05,1.05)

record(f, "X_evolution.gif", profile_names[1:end]; framerate=4) do profile_name
    pname[] = profile_name
end

# ![Movie polytrope](./X_evolution.gif)

##
#=
#### Plot a funny HR diagram

Finally, we can also access the history data of the simulation. We use this to plot a simple HR diagram. As our
microphysics are very simplistic, and the initial condition is not very physical, this looks a bit funny!
=#
f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(T_\text{eff}/[K])", ylabel=L"\log_{10}(L/L_\odot)", xreversed=true)
history = StellarModels.get_history_dataframe_from_hdf5("history.hdf5")
lines!(ax, log10.(history[!, "T_surf"]), log10.(history[!, "L_surf"]))
f

##
#=
#### Perform some cleanup

Internally we want to prevent storing any of the hdf5 files into our git repos, so I remove them. You can also take
advantage of `julia` as a scripting language to post-process your simulation output in a similar way.
=#
rm("history.hdf5")
rm("profiles.hdf5")
rm("example_options.toml")
