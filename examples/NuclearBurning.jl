#=
# NuclearBurning.jl

This notebook provides a simple example of a star with simplified mycrophysics undergoing nuclear burning.
Import all necessary Jems modules. We will also do some benchmarks, so we import BenchmarkTools as well.
=#
using BenchmarkTools
using Jems.Chem
using Jems.Constants
using Jems.EOS
using Jems.Opacity
using Jems.Evolution

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
nvars = 6
nspecies = 2
varnames = [:lnP, :lnT, :lnr, :lum, :H1, :He4]
structure_equations = [Evolution.equationHSE, Evolution.equationT, Evolution.equationContinuity,
                       Evolution.equationLuminosity, Evolution.equationH1, Evolution.equationHe4]
nz = 1000
eos = EOS.IdealEOS(false)
opacity = Opacity.SimpleElectronScatteringOpacity()
sm = StellarModel(varnames, structure_equations, nvars, nspecies, nz, eos, opacity);

##
#=
### Initialize StellarModel and evaluate equations and jacobian

We do not have a working initial condition yet. We require pressure, temperature profiles. One simple available initial
condition is that of an n=1 polytrope. This sets the pressure and density and computes the temperature from the EOS. The
luminosity is initialized by assuming pure radiative transport for the temperature gradient produced by the polytrope.

The normal evolution loop will store the information at the end of the step into an attribute of type `StellarStepInfo`,
stored at `sm.esi`. After initializing our polytrope we can mimic that behavior by calling `set_end_step_info!(sm)`.
TODO: more explanation.
=#
Evolution.n1_polytrope_initial_condition(sm, MSUN, 100 * RSUN; initial_dt=10 * SECYEAR)

Evolution.set_end_step_info!(sm)
Evolution.cycle_step_info!(sm)
Evolution.set_start_step_info!(sm)
Evolution.eval_jacobian_eqs!(sm)

##
#=
### Benchmarking

The previous code leaves everything ready to solve the linearized system. We use the package LinearSolve for this, which
provides various algorithms for linear systems defined by sparse arrays. In the future we might want to try additional
solvers provided (for instance, solvers that make use of the Krylov space such as GMRES). Alternate solvers might work
much better for systems with excessively large nuclear networks.

We can compute a simple benchmark of the time it takes to solve the linear system once. Each timestep will require
multiple iterations of the Newton solver, so this would be a lower bound on the time that will take.
=#
using LinearSolve
@benchmark begin
    $sm.linear_solver.A = $sm.jacobian
    $sm.linear_solver.b = $sm.eqs_numbers
    corr = solve($sm.linear_solver)
end
#=
On my system, this takes on the order of 800 μs. Next up we can check how long it takes to compute a single row of the
jacobian. With a row here I mean all the entries that correspond to one cell. The code below benchmarks the time it
takes to compute the jacobian elements associated with row 2
=#

##

# Benchmark one jacobian row
@benchmark Evolution.eval_jacobian_eqs_row!(sm, 2)

#=
Again on my machine, this takes $\sim 16\;\mathrm{\mu s}$. This is a short amount of time, but we have a thousand cells
to compute. Let's benchmark the calculation of the full jacobian.
=#

##

# Benchmark entire jacobian
@benchmark Evolution.eval_jacobian_eqs!(sm)

#=
And on my computer, this took about $5.2\;\mathrm{ms}$. Even though we have a thousand cells, the computation time was
not a thousand times longer than computing the components of the jacobian for a single cell. The reason for this is that
the calculation is parallelized so cells are done independently. However, I used 8 cores for my calculations, so the
scaling is less than ideal. One of the main culprits here is the garbage collector. Current versions of julia can only
perform garbage collection in a serial way, so it does not take advantage of all threads. Starting with julia 1.10, the
garbage collector will be able to run in multiple threads, so that should alleviate issues with performance scaling.
=#

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
          [solver]
          newton_max_iter_first_step = 1000
          newton_max_iter = 200

          [timestep]
          dt_max_increase = 2.0

          [termination]
          max_model_number = 3000
          max_center_T = 4e7

          [io]
          profile_interval = 50
          """)
end
Evolution.set_options!(sm.opt, "./example_options.toml")
rm(sm.opt.io.hdf5_history_filename; force=true)
rm(sm.opt.io.hdf5_profile_filename; force=true)
Evolution.n1_polytrope_initial_condition(sm, MSUN, 100 * RSUN; initial_dt=1000 * SECYEAR)

@time Evolution.do_evolution_loop(sm)

##
#=
### Plotting with Makie

Now that our simulation is complete we can analyze the results. We make use of the Makie package for this. I'm not a fan
of the Makie defaults, so I adjust them. I normally also adjust the fonts to be consistent with \LaTeX, but I avoid that
here so we don't need to distribute those fonts together with Jems.
=#
using CairoMakie, LaTeXStrings
basic_theme = Theme(
                    #fonts = (regular = "ComputerModernFont/cmunrm.ttf", bold = "ComputerModernFont/cmunbx.ttf", italic = "ComputerModernFont/cmunti.ttf", bold_italic = "ComputerModernFont/cmunbi.ttf"), # taken from https://sourceforge.net/projects/cm-unicode/
                    fontsize=30, resolution=(1000, 750), linewidth=7,
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
profile_names = Evolution.get_profile_names_from_hdf5("profiles.hdf5")

f = Figure()
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(\rho/\mathrm{[g\;cm^{-3}]})", ylabel=L"\log_{10}(P/\mathrm{[dyn]})")

pname = Observable(profile_names[1])

profile = @lift(Evolution.get_profile_dataframe_from_hdf5("profiles.hdf5", $pname))
log10_ρ = @lift($profile[!, "log10_ρ"])
log10_P = @lift($profile[!, "log10_P"])

profile_line = lines!(ax, log10_ρ, log10_P; label="real profile")
xvals = LinRange(-13, 4, 100)
lines!(ax, xvals, (1 + 1 / 1) .* xvals .+ 20; label="n=1")
lines!(ax, xvals, (1 + 1 / (1.5)) .* xvals .+ 15; label="n=1.5")
lines!(ax, xvals, (1 + 1 / 3) .* xvals .+ 15; label="n=3")
axislegend(ax; position=:rb)

model_number_str = @lift("model number=$(parse(Int,$pname))")
profile_text = text!(ax, -10, 20; text=model_number_str)

record(f, "rho_P_evolution.gif", profile_names[1:end]; framerate=2) do profile_name
    pname[] = profile_name
end

# ![Movie polytrope](./rho_P_evolution.gif)

##
#=
### Check nuclear burning

We see that the structure evolves towards an n=3 polytrope. Deviations near the core are due to the non-homogeneous
composition as hydrogen is burnt. We can similarly visualize how the hydrogen mass fraction changes in the simulation.
In here only one frame shows the hydrogen that was burnt, to better visualize that you can adjust `profile_interval` in
the [Io](Evolution.md##Io.jl) options (and probably adjust the framerate).
=#
profile_names = Evolution.get_profile_names_from_hdf5("profiles.hdf5")

f = Figure()
ax = Axis(f[1, 1]; xlabel=L"\mathrm{Mass}\;[M_\odot]", ylabel=L"X")

pname = Observable(profile_names[1])

profile = @lift(Evolution.get_profile_dataframe_from_hdf5("profiles.hdf5", $pname))
mass = @lift($profile[!, "mass"])
X = @lift($profile[!, "X"])
model_number_str = @lift("model number=$(parse(Int,$pname))")

profile_line = lines!(ax, mass, X; label="real profile")
profile_text = text!(ax, 0.7, 0.0; text=model_number_str)

record(f, "X_evolution.gif", profile_names[1:end]; framerate=2) do profile_name
    pname[] = profile_name
end

# ![Movie polytrope](./X_evolution.gif)

##
#=
### Plot a funny HR diagram

Finally, we can also access the history data of the simulation. We use this to plot a simple HR diagram. As our
microphysics are very simplistic, and the initial condition is not very physical, this looks a bit funny!
=#
f = Figure()
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(T_\mathrm{eff}/[K])", ylabel=L"\log_{10}(L/L_\odot)", xreversed=true)
history = Evolution.get_history_dataframe_from_hdf5("history.hdf5")
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
