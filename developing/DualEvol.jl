#=
# NuclearBurning.jl

This notebook provides a simple example of a star with simplified microphysics undergoing nuclear burning.
Import all necessary Jems modules. We will also do some benchmarks, so we import BenchmarkTools as well.
=#
using BenchmarkTools
using ForwardDiff
using Jems.Chem
using Jems.Constants
using Jems.EOS
using Jems.Opacity
using Jems.NuclearNetworks
using Jems.Turbulence
using Jems.StellarModels
using Jems.Evolution
using Jems.Plotting
using Jems.DualSupport

##

net = NuclearNetwork([:H1, :He4], [(:kipp_rates, :kipp_pp)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)

external_tag = ForwardDiff.Tag{:external, nothing}
ForwardDiff.tagcount(external_tag)  # register tag
internal_tag = ForwardDiff.Tag{:internal, nothing}
ForwardDiff.tagcount(internal_tag)
d_to_mass_dual_type = ForwardDiff.Dual{external_tag, Float64, 1}
# d_to_mass_dual_type = Float64

sm = StellarModel(StellarModels.DefaultStellarEquationSet(), nz, nextra, net, eos, opacity, turbulence; number_type=d_to_mass_dual_type, internal_dual_tag=internal_tag);

##

n = 3
mass_dual = ForwardDiff.Dual{external_tag}(1.0 * MSUN, 1.0)
# mass = 5.0 * MSUN
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], mass_dual,
                                             100 * RSUN; initial_dt=10 * SECYEAR)



##
Evolution.compute_starting_model_properties!(sm)

## test some of the dual operations
typeof(sm.props.κ[1])
val00 = get_mixed_00_dual(sm.props.κ[1])
valp1 = get_mixed_p1_dual(sm.props.κ[2])
dm_00 = sm.props.dm[1]
dm_p1 = sm.props.dm[2]
valface_dual = exp((log(val00) * dm_p1 + log(valp1) * dm_00) / (dm_00 + dm_p1))
StellarModels.update_mixed_dual_data!(sm.props.κ_face[1], valface_dual)
StellarModels.eval_mixed_property_log!(sm.props.κ[1], sm.props.κ[2], sm.props.dm[1], sm.props.dm[2], sm.props.κ_face[1])

##
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


##
#set initial condition and run model
n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09],
                                             mass_dual, 100 * RSUN; initial_dt=10 * SECYEAR)
@time Evolution.do_evolution_loop!(sm, plotter=Plotting.NullPlotter());

##
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
ylims!(ax, -0.05, 1.05)

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
