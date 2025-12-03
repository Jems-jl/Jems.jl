#=
# NuclearBurning.jl

This notebook provides a simple example of a star with simplified microphysics undergoing nuclear burning.
Import all necessary Jems modules. We will also do some benchmarks, so we import BenchmarkTools as well.
=#
# using BenchmarkTools
using Jems.Chem
using Jems.Constants
using Jems.EOS
using Jems.Opacity
using Jems.NuclearNetworks 
using Jems.Turbulence
using Jems.StellarModels
using Jems.Evolution
using Jems.Plotting
using Jems.ReactionRates
using Jems.DualSupport


function split_omega(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
    # use same omega in both cells to preserve energy
     varnew_low[sm.vari[:gamma_turb]] = var_00[sm.vari[:gamma_turb]]
     varnew_up[sm.vari[:gamma_turb]] = var_00[sm.vari[:gamma_turb]]
end


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

varnames = [:lnρ, :lnT, :lnr, :lum, :gamma_turb] 
varscaling = [:log, :log, :log, :maxval, :log]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity,
                       Evolution.gammaTurb]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa, split_omega]
net = NuclearNetwork([:H1, :He4], [(:kipp_rates, :kipp_pp)])
# net2 = NuclearNetwork([:H1, :He4, :C12, :N14, :O16], [(:kipp_rates, :kipp_cno)])
# net = merge_nuclear_networks([net1, net2])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.cgMLT(1.0, 9/4)
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
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], 10*MSUN,
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
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
Evolution.eval_jacobian_eqs!(sm)
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
          newton_max_iter = 200
          scale_max_correction = 1.0
          solver_progress_iter = 1
          relative_correction_tolerance = 1e12

          [timestep]
          dt_max_increase = 1.5
          delta_R_limit = 0.01
          delta_Tc_limit = 0.01     
          delta_Xc_limit = 0.005

          [termination]
          max_model_number = 100000
          max_center_T = 1e8
          

          
          [io]
          profile_interval = 50
          terminal_header_interval = 100
          terminal_info_interval = 100
          profile_values = ["zone", "mass", "dm", "log10_rho", "log10_r", "log10_P", "log10_T", "luminosity",
                                      "X", "Y","D_face", "nabla_a_face", "nabla_r_face","nabla_face", "nabla_tdc","velocity_turb","turb_energy", "D_face_kuhfuss"]
            history_values = ["age", "dt", "star_mass", "alpha_overshoot"] 
          """)
end
StellarModels.set_options!(sm.opt, "./example_options.toml")
rm(sm.opt.io.hdf5_history_filename; force=true)
rm(sm.opt.io.hdf5_profile_filename; force=true)

#Configure live plots. To turn off one can use `plotter = Plotting.NullPlotter()`
using GLMakie
set_theme!(Plotting.basic_theme())
f = Figure(size=(1400,750))
hist_plot = Plotting.HistoryPlot(f[1,3], sm, x_name="age", y_name="alpha_overshoot", link_yaxes=true)
ylims!(hist_plot.axis, 0, 0.5)
plots = [Plotting.HRPlot(f[1,1]),
         Plotting.TRhoProfile(f[1,2]),
         Plotting.KippenLine(f[2,1], xaxis=:time, time_units=:Gyr),
         Plotting.AbundancePlot(f[2,2],net,log_yscale=true, ymin=1e-3),
         Plotting.HistoryPlot(f[3,1], sm, x_name="age", y_name="X_center", othery_name="Y_center", link_yaxes=true),
         hist_plot,
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


# Output directory
outdir = "/home/ritavash/Desktop/Resources/convection_results/"

# Loop over models 100 to 900 in steps of 50
for model_number in 100:50:500
    # Zero-padded model name
    model_str = lpad(model_number, 10, '0')

    # Load profile
    profile = StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", model_str)

    # Select core region (mass ≤ 1 M_sun)
    core_profile = profile[profile[!, "mass"] .<= 1, :]

    # Create figure
    f = Figure(resolution = (1200, 600))
    ax = Axis(f[1, 1];
        xlabel = L"Mass\,[M_\odot]",
        ylabel = L"v_\mathrm{turb}",
        title = "Turbulent velocity — Model $(model_number)"
    )

    # Plot turbulent velocity
    lines!(ax,
        core_profile[!, "mass"],
        core_profile[!, "velocity_turb"];
        color = :blue,
        linewidth = 3
    )

    # Output filename
    outfile = outdir * "kuhfuss_tdd_" * string(model_number) * ".png"

    # Save figure
    save(outfile, f)
    println("Saved: $outfile")
end
##

#  save("/home/ritavash/Desktop/Resources/convection_results/D_kuhfuss_vs_D_mlt500.png", f)

##
# Plotting v_turb_kuhfuss vs v_steady

f= Figure(resolution = (1200, 800));
profile = StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", "0000000500")
core_profile = profile[profile[!, "mass"] .<= 100, :]
ax1 = Axis(f[1,1];xlabel = L"Mass\;[M_\odot]", ylabel = "Velocity gradients", title = "profiles")

lines!(ax1,
    core_profile[!, "mass"],
    ((core_profile[!, ("velocity_turb")])) ;
    label = L"v_{\mathrm{turb,\,Kuhfuss}}", color = :red, linewidth = 3
)

lines!(ax1,
    core_profile[!, "mass"],
    core_profile[!, "v_face"] ;
    label = L"v_{mlt}", color = :green, linewidth = 3
)



axislegend(ax1, position = :rt) 
f
save("/home/ritavash/Desktop/Resources/convection_results/velocity_gradients.png", f)
##
f= Figure(resolution = (1200, 800));
profile = StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", "0000000400")



core_profile = profile[profile[!, "mass"] .<= 0.4, :]
ax1 = Axis(f[1,1];xlabel = L"Mass\;[M_\odot]", ylabel = "Gradients", title = "D_mixing at the core of the star")
lines!(ax1,
    core_profile[!, "mass"],
    core_profile[!, "D_face"];
    label = L"D_\mathrm{turb}", color = :blue, linewidth = 3
)
##
f= Figure(resolution = (1400, 1000));
profile = StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", "0000000400")
core_profile = profile[0.207 .<=profile[!, "mass"] .<= 1, :]
ax1 = Axis(f[1,1];xlabel = L"Mass\;[M_\odot]", ylabel = "Gradients", title = "Temperature Gradient Comparison at the core of the star")
lines!(ax1,
    core_profile[!, "mass"],
    core_profile[!, "nabla_face"];
    label = L"\nabla_\mathrm{actual}", color = :blue, linewidth = 3
)
lines!(ax1,
    core_profile[!, "mass"],
    core_profile[!, "nabla_r_face"];
    label = L"\nabla_\mathrm{rad}", color = :red, linewidth = 3, linestyle = :dash
)
lines!(ax1,
    core_profile[!, "mass"],
    core_profile[!, "nabla_a_face"];
    label = L"\nabla_\mathrm{ad}", color = :green, linewidth = 3, linestyle = :dot
)           
axislegend(ax1; position = :rt)
f
# save("/home/ritavash/Desktop/Resources/convection_results/Temperature_gradient_kuhfuss.png", f)
##

# 1. Setup and Data Filtering
f = Figure(resolution = (1400, 1500));
profile = StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", "0000000600")

"""
Calculating the overshooting length 
"""
m_profile = profile[!, "mass"]
P_profile = profile[!, "log10_P"]
ρ_profile = profile[!, "log10_rho"]
r_profile = profile[!, "log10_r"]
nabla_rad = profile[!, "nabla_r_face"]
nabla_ad = profile[!, "nabla_a_face"]
#get the start of the convective boundary/ Schwarzschild boundary.
boundary_indices = findall(nabla_ad .> nabla_rad)

function interpolate_sch_radius(profile, boundary_index)

    i_rad = boundary_index
    i_conv = boundary_index 

    X_rad = profile[i_rad, "log10_r"]
    X_conv = profile[i_conv, "log10_r"]

    dnabla_rad  = profile[i_rad, "nabla_r_face"] - profile[i_rad, "nabla_a_face"]
    dnabla_conv = profile[i_conv, "nabla_r_face"] - profile[i_conv, "nabla_a_face"]

    dnabla_total = dnabla_rad - dnabla_conv

    if dnabla_total == 0.0
        # Should not happen unless profiles are flat, return the midpoint
        logr_sch = (X_rad + X_conv) / 2.0
    else
        # f = -Δ∇_conv / (Δ∇_rad - Δ∇_conv)
        logr_sch = X_conv - dnabla_conv * (X_rad - X_conv) / dnabla_total
    end
    return logr_sch
end

if !isempty(boundary_indices)
    boundary_index  = boundary_indices[1]
    m_boundary = m_profile[boundary_index]*MSUN
    P_boundary = 10^(P_profile[boundary_index])
    ρ_boundary = exp(ρ_profile[boundary_index]/log10(ℯ))
    r = interpolate_sch_radius(profile, boundary_index) # r is the interpolated sch boundary 
    r_boundary = exp(( r + log10(RSUN))/log10(ℯ))
    scale_height = P_boundary /(ρ_boundary * CGRAV * m_boundary / r_boundary^2)
    scaled_radius = r_boundary / RSUN 
end    
scale_height_arr = []
for i in eachindex(m_profile)
    radius = exp((r_profile[i]+ log10(RSUN))/log10(ℯ))
    scaled_measure = radius - scaled_radius
    push!(scale_height_arr, scaled_measure)
end

radiative_zone_indices = findall((profile[!, "D_face_kuhfuss"]) .< 1e5)

if !isempty(radiative_zone_indices)
    radiative_zone_index = radiative_zone_indices[1]
    radiative_boundary_mass = profile[radiative_zone_index, "mass"]
    r_overshoot = 10.0^(r_profile[radiative_zone_index]+ log10(RSUN))
    overshooting_length = r_overshoot - r_boundary
    return overshooting_length / scale_height

end
alpha_overshoot = overshooting_length / scale_height
core_profile = profile[0 .<=profile[!, "mass"] .<= 0.5, :]

# 2. Find the Convective Boundary
nabla_rad = core_profile[!, "nabla_r_face"]
nabla_ad = core_profile[!, "nabla_a_face"]

# Find the first mass point (starting from the center) where the convection condition (nabla_rad > nabla_ad) fails.
boundary_indices = findall(nabla_ad .>= nabla_rad)

boundary_mass = nothing
if !isempty(boundary_indices)
    # Use the mass coordinate at the boundary index
    boundary_index = boundary_indices[1]
    boundary_mass = core_profile[boundary_index, "mass"]
end

# 3. Plotting
ax1 = Axis(f[1,1];
    xlabel = L"Mass\;[M_\odot]",
    ylabel = " log v_turb",
    title = "v_turb at the Core of the Star (1M☉)"
)
ax3 = Axis(f[2,1];
    xlabel = L"Mass\;[M_\odot]",
    ylabel = " log D_mixing",
    title = "D_mixing at the Core of the Star (1M☉)",
)
ax2 = Axis(f[3,1];
    xlabel = L"Mass\;[M_\odot]",
    ylabel = "Gradients",
    title = "Gradients at the Core of the Star (1M☉)"
)
# ylims!(ax2,0.3993542, 0.3995)
# xlims!(ax2,0.23, 0.37)
# Plot D_mixing (on a log scale)
lines!(ax1,
    core_profile[!, "mass"],
    log10.(core_profile[!, "velocity_turb"]) ;
    label = L"log_{10}(v_{\mathrm{turb,\,Kuhfuss}})", color = :red, linewidth = 3
)
lines!(ax3,
    core_profile[!, "mass"],
    log10.(core_profile[!, "D_face_kuhfuss"]) ;
    label = L"log_{10}(D_{\mathrm{turb,\,Kuhfuss}})", color = :blue, linewidth = 3
)


# Plot Radiative Gradient
lines!(ax2,
    core_profile[!, "mass"],
    nabla_rad;
    label = L"\nabla_\mathrm{rad}", color = :blue, linewidth = 3, linestyle = :dash
)

# Plot Adiabatic Gradient
lines!(ax2,
    core_profile[!, "mass"],
    nabla_ad;
    label = L"\nabla_\mathrm{ad}", color = :green, linewidth = 3, linestyle = :dot
)
lines!(ax2,
    core_profile[!, "mass"],
    core_profile[!, "nabla_tdc"];
    label = L"\nabla_\mathrm{actual}", color = :red, linewidth = 3
)
# 4. Add the Vertical Line
if boundary_mass !== nothing
    vlines!(ax1, [boundary_mass];
        color = :black,
        linestyle = :dash,
        linewidth = 0.5,
        label = "Convective Boundary"
    )
end
if boundary_mass !== nothing
    vlines!(ax2, [boundary_mass];
        color = :black,
        linestyle = :dash,
        linewidth = 0.5,
        label = "Convective Boundary"
    )
end
if boundary_mass !== nothing
    vlines!(ax3, [boundary_mass];
        color = :black,
        linestyle = :dash,
        linewidth = 0.5,
        label = "Convective Boundary"
    )
end
if boundary_mass !== nothing
    vlines!(ax3, [radiative_boundary_mass];
        color = :black,
        linestyle = :dash,
        linewidth = 0.5,
        label = "Overshooting Boundary"
    )
end
if !isnan(alpha_overshoot)
    # Add alpha_overshoot value to the D_mixing plot (ax3)
    # Position: top-right corner, converting the alpha value to a formatted string.
    text!(ax3, 
        0.31, # X position (mass)
        log10(100), # Y position (log10(D_turb)) - placing it near the top of the plot
        text = L"\alpha_{\mathrm{ov}} = %$(round(alpha_overshoot, digits=3)) H_P",
        align = (:right, :center),
        fontsize = 20,
        color = :black
    )
end

axislegend(ax1, position = :rt)
axislegend(ax2, position = :lt,labelsize=18)

f
# save("/home/ritavash/Desktop/Resources/convection_results/turb_plots_kuhfuss_overshoot_length_1M.png", f)







##
f= Figure(resolution = (1200, 800));
profile = StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", "0000000500")  
core_profile = profile[profile[!, "mass"] .<= 0.1, :]
ax1 = Axis(f[1,1];xlabel = L"Mass\;[M_\odot]", ylabel = "Convective velocity", title = "profiles (model 400)", yscale = log10)

lines!(ax1,
    core_profile[!, "mass"],
    ((core_profile[!, ("velocity_turb")])) ;
    label = L"v_{\mathrm{turb,\,Kuhfuss}}", color = :red, linewidth = 3
)

f
# save("/home/ritavash/Desktop/Resources/convection_results/v_turb_kuhfuss.png", f)

##
profile = StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", "0000000500")

m_profile = profile[!, "mass"]
P_profile = profile[!, "log10_P"]
ρ_profile = profile[!, "log10_ρ"]
r_profile = profile[!, "log10_r"]
nabla_rad = profile[!, "nabla_r_face"]
nabla_ad = profile[!, "nabla_a_face"]
#get the start of the convective boundary/ Schwarzschild boundary.
boundary_indices = findall(nabla_ad .> nabla_rad)

if !isempty(boundary_indices)
    boundary_index  = boundary_indices[1]
    m_boundary = m_profile[boundary_index]*MSUN
    P_boundary = P_profile[boundary_index]/log10(ℯ)
    ρ_boundary = exp(ρ_profile[boundary_index]/log10(ℯ))
    r_boundary = exp((r_profile[boundary_index]+ log10(RSUN))/log10(ℯ))
    scale_height = P_boundary /(ρ_boundary * CGRAV * m_boundary / r_boundary^2)
end  


radiative_zone_indices = findall((profile[!, "D_face_kuhfuss"]) .< 1e-1)

if !isempty(radiative_zone_indices)
    radiative_zone_index = radiative_zone_indices[1]
    r_overshoot = 10.0^(r_profile[radiative_zone_index]+ log10(RSUN))
    overshooting_length = r_overshoot - r_boundary
    return overshooting_length / scale_height

end
alpha_overshoot = overshooting_length / scale_height
println("The overshooting parameter is: ", alpha_overshoot)
println("The overshooting length index is: ", radiative_zone_index)
println("The sch boundary index is: ", boundary_index)
##

"""
Visualizing the overshooting length in a different way 
"""

f = Figure(resolution = (1400, 900));
profile = StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", "0000000600")

"""
Calculating the overshooting length 
"""
m_profile = profile[!, "mass"]
P_profile = profile[!, "log10_P"]
ρ_profile = profile[!, "log10_rho"]
r_profile = profile[!, "log10_r"]
nabla_rad = profile[!, "nabla_r_face"]
nabla_ad = profile[!, "nabla_a_face"]
#get the start of the convective boundary/ Schwarzschild boundary.
boundary_indices = findall(nabla_ad .> nabla_rad)

function interpolate_sch_radius(profile, boundary_index)

    i_rad = boundary_index
    i_conv = boundary_index 

    X_rad = profile[i_rad, "log10_r"]
    X_conv = profile[i_conv, "log10_r"]

    dnabla_rad  = profile[i_rad, "nabla_r_face"] - profile[i_rad, "nabla_a_face"]
    dnabla_conv = profile[i_conv, "nabla_r_face"] - profile[i_conv, "nabla_a_face"]

    dnabla_total = dnabla_rad - dnabla_conv

    if dnabla_total == 0.0
        # Should not happen unless profiles are flat, return the midpoint
        logr_sch = (X_rad + X_conv) / 2.0
    else
        # f = -Δ∇_conv / (Δ∇_rad - Δ∇_conv)
        logr_sch = X_conv - dnabla_conv * (X_rad - X_conv) / dnabla_total
    end
    return logr_sch
end

if !isempty(boundary_indices)
    boundary_index  = boundary_indices[1]
    m_boundary = m_profile[boundary_index]*MSUN
    P_boundary = 10^(P_profile[boundary_index])
    ρ_boundary = exp(ρ_profile[boundary_index]/log10(ℯ))
    r = interpolate_sch_radius(profile, boundary_index) # r is the interpolated sch boundary 
    r_boundary = exp(( r + log10(RSUN))/log10(ℯ))
    scale_height = P_boundary /(ρ_boundary * CGRAV * m_boundary / r_boundary^2)
    scaled_radius = r_boundary / RSUN 
end   
scale_height_arr = []
for i in eachindex(m_profile)
    radius = exp((r_profile[i]+ log10(RSUN))/log10(ℯ))
    scaled_measure = radius - scaled_radius
    push!(scale_height_arr, scaled_measure)
end


radiative_zone_indices = findall((profile[!, "D_face_kuhfuss"]) .< 1e5)

if !isempty(radiative_zone_indices)
    radiative_zone_index = radiative_zone_indices[1]
    radiative_boundary_mass = profile[radiative_zone_index, "mass"]
    r_overshoot = 10.0^(r_profile[radiative_zone_index]+ log10(RSUN))
    overshooting_length = r_overshoot - r_boundary
    return overshooting_length / scale_height

end
alpha_overshoot = overshooting_length / scale_height
core_profile = profile[0 .<=profile[!, "mass"] .<= 0.5, :]

# 2. Find the Convective Boundary
nabla_rad = core_profile[!, "nabla_r_face"]
nabla_ad = core_profile[!, "nabla_a_face"]


ax3 = Axis(f[1,1];
    xlabel = L"r(m) - r_{\mathrm{sch}}/H_{P,sch}",
    ylabel = "D_mix",
    title = "D mixing at the Core of the Star: M500 (1M☉)",
)

lines!(ax3,
    scale_height_arr[0 .<=profile[!, "mass"] .<= 0.5],
    log10.(core_profile[!, "D_face_kuhfuss"]) ;
    label = L"log_{10}(D_{\mathrm{turb,\,Kuhfuss}})", color = :blue, linewidth = 3
)
# if boundary_mass !== nothing
#     vlines!(ax3, scale_height_arr[boundary_index];
#         color = :black,
#         linestyle = :dash,
#         linewidth = 0.5,
#         label = "Convective Boundary"
#     )
# end
# if boundary_mass !== nothing
#     vlines!(ax3, scale_height_arr[radiative_zone_index];
#         color = :black,
#         linestyle = :dash,
#         linewidth = 0.5,
#         label = "Overshooting Boundary"
#     )
# end
if !isnan(alpha_overshoot)
    text!(ax3, 
        scale_height_arr[boundary_index], # X position (mass)
        log10(100), # Y position (log10(D_turb)) - placing it near the top of the plot
        text = L"\alpha_{\mathrm{ov}} = %$(round(alpha_overshoot, digits=3)) H_P",
        align = (:right, :center),
        fontsize = 20,
        color = :black
    )
end
f
# save("/home/ritavash/Desktop/Resources/convection_results/D_mix_visualization_1M.png", f)
##
"""
Plot for Visualizing the temperature gradients 
"""

f= Figure(resolution = (1400, 1000));
profile = StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", "0000000500")
"""
Calculating the overshooting length 
"""
m_profile = profile[!, "mass"]
P_profile = profile[!, "log10_P"]
ρ_profile = profile[!, "log10_ρ"]
r_profile = profile[!, "log10_r"]
nabla_rad = profile[!, "nabla_r_face"]
nabla_ad = profile[!, "nabla_a_face"]
#get the start of the convective boundary/ Schwarzschild boundary.
boundary_indices = findall(nabla_ad .>= nabla_rad)
function interpolate_sch_radius(profile, boundary_index)

    i_rad = boundary_index
    i_conv = boundary_index - 1

    X_rad = profile[i_rad, "log10_r"]
    X_conv = profile[i_conv, "log10_r"]

    dnabla_rad  = profile[i_rad, "nabla_r_face"] - profile[i_rad, "nabla_a_face"]
    dnabla_conv = profile[i_conv, "nabla_r_face"] - profile[i_conv, "nabla_a_face"]

    dnabla_total = dnabla_rad - dnabla_conv

    if dnabla_total == 0.0
        # Should not happen unless profiles are flat, return the midpoint
        logr_sch = (X_rad + X_conv) / 2.0
    else
        # f = -Δ∇_conv / (Δ∇_rad - Δ∇_conv)
        logr_sch = X_conv - dnabla_conv * (X_rad - X_conv) / dnabla_total
    end
    return logr_sch
end

if !isempty(boundary_indices)
    boundary_index  = boundary_indices[1]
    m_boundary = m_profile[boundary_index]*MSUN
    boundary_mass = m_profile[boundary_index]
    P_boundary = P_profile[boundary_index]/log10(ℯ)
    ρ_boundary = exp(ρ_profile[boundary_index]/log10(ℯ))
    r = interpolate_sch_radius(profile, boundary_index) # r is the interpolated sch boundary 
    r_boundary = exp(( r + log10(RSUN))/log10(ℯ))
    scale_height = P_boundary /(ρ_boundary * CGRAV * m_boundary / r_boundary^2)
    scaled_radius = r_boundary / RSUN 
end   

core_profile = profile[0.1.<=profile[!, "mass"] .<= 0.5, :]
nabla_diff = core_profile[!, "nabla_tdc"] - core_profile[!, "nabla_a_face"]
ax1 = Axis(f[1,1];xlabel = L"Mass\;[M_\odot]", ylabel = "Gradients", title = "Temperature Gradient Comparison at the core of the star")
ylims!(ax1,-0.2, 0.2)
lines!(ax1,
    core_profile[!, "mass"],
    nabla_diff;
    label = L"\nabla - \nabla_\mathrm{ad}", color = :blue, linewidth = 3
)
if boundary_mass !== nothing
    vlines!(ax1, [boundary_mass];
        color = :black,
        linestyle = :dash,
        linewidth = 0.8,
        label = "Convective Boundary"
    )
end
axislegend(ax1; position = :rt)
f
save("/home/ritavash/Desktop/Resources/convection_results/temperature_gradient_difference_1.png", f)



##
history = StellarModels.get_history_dataframe_from_hdf5("history.hdf5")

# 2. Create the figure and axis
f = Figure(resolution = (1400, 1000));
ax = Axis(f[1, 1]; 
          xlabel=L"\mathrm{Age}\;[\mathrm{yr}]", 
          ylabel=L"\alpha_{\mathrm{ov}}/\mathrm{H_p}", 
          title=L"\text{Convective Overshoot Parameter History (1} M_\odot)"
)
# lines!(ax, 
#     history[!, "age"], 
#     history[!, "alpha_overshoot"], 
#     # Line-specific parameters:
#     linewidth = 3, 
#     color = :navy, 
#     label = L"\alpha_{\mathrm{ov}}"
# )
# 3. Plot the history (Age vs. Alpha_Overshoot)
# scatter!(ax, 
#     history[!, "age"], 
#     history[!, "alpha_overshoot"], 
#     # Scatter-specific parameters:
#     markersize = 8, 
#     color = :firebrick, 
#     label = L"\alpha_{\mathrm{ov}}"
# )
scatterlines!(ax, history[!, "age"], 
    history[!, "alpha_overshoot"],
    linewidth = 0.5)
 ylims!(ax, 0, 0.2)
f
save("/home/ritavash/Desktop/Resources/convection_results/alpha_overshoot/1_history.png", f)
##
f= Figure(resolution = (1200, 800));
profile = StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", "0000000600")
core_profile = profile[profile[!, "mass"] .<= 1, :]
ax1 = Axis(f[1,1];xlabel = L"Mass\;[M_\odot]", ylabel = "Turbulent energy", title = "profiles")

scatterlines!(ax1,
    core_profile[!, "mass"],
    core_profile[!, "turb_energy"];
    label = L"E_{\mathrm{turb,\,Kuhfuss}}", color = :red, linewidth = 0.5
)



axislegend(ax1, position = :rt) 
f
# save("/home/ritavash/Desktop/Resources/convection_results/velocity_gradients.png", f)


##
#  history[!,"i_ov"][500]
##
### Perform some cleanup
#=
Internally we want to prevent storing any of the hdf5 files into our git repos, so I remove them. You can also take
advantage of `julia` as a scripting language to post-process your simulation output in a similar way.
=#
rm("history.hdf5")
rm("profiles.hdf5")
rm("example_options.toml")



