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
using Jems.ReactionRates
using Jems.DualSupport
using CairoMakie

function split_omega(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
    # use same omega in both cells to preserve energy
     varnew_low[sm.vari[:gamma_turb]] = var_00[sm.vari[:gamma_turb]]
     varnew_up[sm.vari[:gamma_turb]] = var_00[sm.vari[:gamma_turb]]
end


function calculate_nabla_tdc(sm:: StellarModel, k :: Int)
    L = get_00_dual(sm.props.L[k]) * LSUN
    γ₀ = get_00_dual(sm.props.gamma_turb[k])
    ω = exp(γ₀)
    m₀ = sm.props.m[k]
    r₀ = exp(get_00_dual(sm.props.lnr[k]))
    if k == sm.props.nz
        P = get_00_dual(sm.props.eos_res[k].P)
        ρ = get_00_dual(sm.props.eos_res[k].ρ)
        T = get_00_dual(sm.props.eos_res[k].T)
        κ = get_00_dual(sm.props.κ[k])
        cₚ = get_00_dual(sm.props.eos_res[k].cₚ)
        Hₚ = P / (ρ * CGRAV * m₀ / r₀^2)
        Λ = 1/(1/Hₚ + 1/r₀)
        k_rad = 16 * SIGMA_SB * T^3 / (3 * κ * ρ)
        α₂ = ρ*cₚ*0.5*sqrt(2/3)*Λ*sqrt(ω)
        ∇ᵣ = 3 * κ * L * P / (16π * CRAD * CLIGHT * CGRAV * m₀ * T^4)
        ∇ₐ = get_00_dual(sm.props.eos_res[k].∇ₐ)
        SA = (∇ᵣ - ∇ₐ)*(1 + α₂/k_rad)^(-1)        
        ∇ = ∇ₐ + SA
        return ∇
    end 

    P_face = exp(get_00_dual(sm.props.lnP_face[k]))
    ρ_face = exp(get_00_dual(sm.props.lnρ_face[k]))
    T_face = exp(get_00_dual(sm.props.lnT_face[k]))
    cₚ =  get_00_dual(sm.props.cₚ_face[k])
    κ = get_00_dual(sm.props.κ_face[k])
    Hₚ = P_face / (ρ_face * CGRAV * m₀ / r₀^2) 
    Λ = 1/(1/Hₚ + 1/r₀)
    m₀ = sm.props.m[k]
    k_rad = 16 * SIGMA_SB * T_face^3 / (3 * κ * ρ_face)
    α₂ = ρ_face*cₚ*0.5*sqrt(2/3)*Λ*sqrt(ω)
    ∇ᵣ = 3 * κ * L * P_face / (16π * CRAD * CLIGHT * CGRAV * m₀ * T_face^4)
    ∇ₐ = get_00_dual(sm.props.∇ₐ_face[k])
    SA = (∇ᵣ - ∇ₐ)*(1 + α₂/k_rad)^(-1)
    ∇ = ∇ₐ + SA 
    return ∇
end 
function get_nabla_tdc(sm, k)
    # Calls the calculation function and extracts the pure numerical value.
    return calculate_nabla_tdc(sm, k).value
end


function interpolate_schwarzschild_boundary_live(sm::StellarModel, i_rad::Int)
    i_conv = i_rad - 1
    
    #Radius 
    logr_rad  = get_value(sm.props.lnr[i_rad])
    logr_conv = get_value(sm.props.lnr[i_conv])
    #Pressure
    logP_rad  = get_value(sm.props.eos_res[i_rad].P)
    logP_conv = get_value(sm.props.eos_res[i_conv].P)
    #density 
    logρ_rad  = get_value(sm.props.lnρ[i_rad])
    logρ_conv = get_value(sm.props.lnρ[i_conv])
    #mass
    m_rad  = sm.props.m[i_rad] 
    m_conv = sm.props.m[i_conv]

    # Delta-Nabla (Δ∇ = ∇_rad - ∇_ad)
    dnabla_rad  = get_value(sm.props.turb_res[i_rad].∇ᵣ) - get_value(sm.props.∇ₐ_face[i_rad])
    dnabla_conv = get_value(sm.props.turb_res[i_conv].∇ᵣ) - get_value(sm.props.∇ₐ_face[i_conv])
    
    # 2. Linear Interpolation Calculation
    dnabla_total = dnabla_rad - dnabla_conv
    
    if abs(dnabla_total) < 1e-12 
        # Avoid division by zero: return the midpoint if profiles are flat
        logr_sch = (logr_rad + logr_conv) / 2.0
        logP_sch = (logP_rad + logP_conv) / 2.0
        logρ_sch = (logρ_rad + logρ_conv) / 2.0
        m_sch = (m_rad + m_conv) / 2.0
    else
        # Interpolation: logr_sch = logr_conv - Δ∇_conv * (Δlogr / Δ(Δ∇))
        logr_sch = logr_conv - dnabla_conv * (logr_rad - logr_conv) / dnabla_total
        logP_sch = logP_conv - dnabla_conv * (logP_rad - logP_conv) / dnabla_total
        logρ_sch = logρ_conv - dnabla_conv * (logρ_rad - logρ_conv) / dnabla_total
        m_sch = m_conv - dnabla_conv * (m_rad - m_conv) / dnabla_total
    end
    
    # Return the interpolated natural log radius (ln r_sch)
    return logr_sch, logP_sch, logρ_sch, m_sch
end

function calculate_overshoot_length(sm:: StellarModel) 
    

    OVERSHOOT_THRESHOLD = 1e-1 
    
    nz_interior  = sm.props.nz - 1 
    i_Sch = nothing
    
    # 1. Find Schwarzschild Boundary Index (i_Sch: first point where ∇_ad > ∇_rad)
    for k in 2:nz_interior
        nabla_ad = get_value(sm.props.∇ₐ_face[k])
        nabla_rad = get_value(sm.props.turb_res[k].∇ᵣ)

        if nabla_ad > nabla_rad
            i_Sch = k
            break
        end   
    end 
    
    if isnothing(i_Sch)
         return 0.0
    end
    
    # 2. Calculate HP_Sch (using properties at the boundary index)
        # Interpolate to find boundary properties
        logr_sch, logP_sch, logρ_sch, m_sch = interpolate_schwarzschild_boundary_live(sm, i_Sch)
    P_boundary = logP_sch
    ρ_boundary = exp(logρ_sch)
    r_boundary = exp(logr_sch) # Radius in cm
    m_boundary = m_sch
    
    HP_Sch = P_boundary / (ρ_boundary * CGRAV * m_boundary / r_boundary^2)

    # 3. Find Overshoot Edge (i_ov: first point where D_turb is negligible)
    i_ov = nothing

    for k = i_Sch + 1:nz_interior
        D_turb = get_value(sm.props.D_turb[k])
        if D_turb < OVERSHOOT_THRESHOLD
            i_ov = k
            break
        end
    end
    
    if isnothing(i_ov)
        # If mixing doesn't fall below threshold before the surface
        return 0.0 
    end
    
    # 4. Calculate Alpha_ov
    r_overshoot = exp(get_value(sm.props.lnr[i_ov])) # Radius in cm
    
    overshooting_distance_cm = abs(r_overshoot - r_boundary)
    alpha_ov = overshooting_distance_cm / HP_Sch
    
    return alpha_ov 
end

function get_overshoot(sm)
    # Calls the calculation function and extracts the pure numerical value.
    return calculate_overshoot_length(sm)
end

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

varnames = [:lnρ, :lnT, :lnr, :lum, :gamma_turb] 
varscaling = [:log, :log, :log, :maxval, :log]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity,
                       Evolution.gammaTurb]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa, split_omega]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

The normal evolution loop will store the information at the end of the step into an attribute of type `StellarStepInfo`,
stored at `sm.esi` (_end step info_). After initializing our polytrope we can mimic that behavior by calling 
`set_end_step_info!(sm)`. We then 'cycle' this info into the information of a hypothetical previous step with
`cycle_step_info`, so now `sm.psi` contains our initial condition. Finally we call `set_start_step_info` to use `sm.psi`
(_previous step info_) to populate the information needed before the Newton solver in `sm.ssi` (_start step info_).
At last we are in position to evaluate the equations and compute the Jacobian.
=#
n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], 10*MSUN,
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
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
Evolution.eval_jacobian_eqs!(sm)
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
          profile_values = ["zone", "mass", "dm", "log10_ρ", "log10_r", "log10_P", "log10_T", "luminosity",
                                      "X", "Y", "velocity_turb", "omega_e","v_face","D_face", "D_face_kuhfuss", "nabla_a_face", "nabla_r_face","nabla_face", "nabla_tdc"]
            history_values = ["age", "dt", "star_mass", "alpha_overshoot"] 
          """)
end
StellarModels.set_options!(sm.opt, "./example_options.toml")
rm(sm.opt.io.hdf5_history_filename; force=true)
rm(sm.opt.io.hdf5_profile_filename; force=true)

add_profile_option("velocity_turb", "unitless", (sm, k) -> sqrt(2*exp(get_value(sm.props.gamma_turb[k]))))
add_profile_option("omega_e", "unitless", (sm, k) -> sqrt(get_value(sm.props.eos_res[k].P)*1e-14/get_value(sm.props.eos_res[k].ρ)))
add_profile_option("v_face", "unitless", (sm, k) -> get_value(sm.props.turb_res[k].v_turb))
add_profile_option("vel", "unitless", (sm, k) -> sqrt(2*exp(get_value(sm.props.gamma_turb[k]))))
# add_profile_option("nabla_kuf", "unitless", (sm, k) -> get_value(sm.props.∇_face[k]))
add_profile_option("nabla_tdc", "unitless", get_nabla_tdc)
add_history_option("alpha_overshoot", "H_p", get_overshoot)
# add_history_option("i_Sch", "index", get_i_Sch)
# add_history_option("i_ov", "index", get_i_ov)
n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], 
                                            16 * MSUN, 30 * RSUN; initial_dt=10 * SECYEAR)
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
boundary_indices = findall(nabla_ad .> nabla_rad)

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
    P_boundary = P_profile[boundary_index]/log10(ℯ)
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

radiative_zone_indices = findall((profile[!, "D_face_kuhfuss"]) .< 1e-1)

if !isempty(radiative_zone_indices)
    radiative_zone_index = radiative_zone_indices[1]
    radiative_boundary_mass = profile[radiative_zone_index, "mass"]
    r_overshoot = 10.0^(r_profile[radiative_zone_index]+ log10(RSUN))
    overshooting_length = r_overshoot - r_boundary
    return overshooting_length / scale_height

end
alpha_overshoot = overshooting_length / scale_height
core_profile = profile[0.1 .<=profile[!, "mass"] .<= 15, :]

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
    title = "v_turb at the Core of the Star (16M☉)"
)
ax3 = Axis(f[2,1];
    xlabel = L"Mass\;[M_\odot]",
    ylabel = " log D_mixing",
    title = "D_mixing at the Core of the Star (16M☉)",
)
ax2 = Axis(f[3,1];
    xlabel = L"Mass\;[M_\odot]",
    ylabel = "Gradients",
    title = "Gradients at the Core of the Star (16M☉)"
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
save("/home/ritavash/Desktop/Resources/convection_results/turb_plots_kuhfuss_overshoot_length_16M.png", f)







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
boundary_indices = findall(nabla_ad .> nabla_rad)

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
    P_boundary = P_profile[boundary_index]/log10(ℯ)
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


radiative_zone_indices = findall((profile[!, "D_face_kuhfuss"]) .< 1e-1)

if !isempty(radiative_zone_indices)
    radiative_zone_index = radiative_zone_indices[1]
    radiative_boundary_mass = profile[radiative_zone_index, "mass"]
    r_overshoot = 10.0^(r_profile[radiative_zone_index]+ log10(RSUN))
    overshooting_length = r_overshoot - r_boundary
    return overshooting_length / scale_height

end
alpha_overshoot = overshooting_length / scale_height
core_profile = profile[0 .<=profile[!, "mass"] .<= 15, :]

# 2. Find the Convective Boundary
nabla_rad = core_profile[!, "nabla_r_face"]
nabla_ad = core_profile[!, "nabla_a_face"]


ax3 = Axis(f[1,1];
    xlabel = L"r(m) - r_{\mathrm{sch}}/H_{P,sch}",
    ylabel = "D_mix",
    title = "D mixing at the Core of the Star: M500 (16M☉)",
)

lines!(ax3,
    scale_height_arr[0 .<=profile[!, "mass"] .<= 15],
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
save("/home/ritavash/Desktop/Resources/convection_results/D_mix_visualization_16M.png", f)
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

core_profile = profile[0.5.<=profile[!, "mass"] .<= 15, :]
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
save("/home/ritavash/Desktop/Resources/convection_results/temperature_gradient_difference_16.png", f)



##
history = StellarModels.get_history_dataframe_from_hdf5("history.hdf5")

# 2. Create the figure and axis
f = Figure(resolution = (1400, 1000));
ax = Axis(f[1, 1]; 
          xlabel=L"\mathrm{Age}\;[\mathrm{yr}]", 
          ylabel=L"\alpha_{\mathrm{ov}}/\mathrm{H_p}", 
          title=L"\text{Convective Overshoot Parameter History (16} M_\odot)"
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
scatter!(ax, 
    history[!, "age"], 
    history[!, "alpha_overshoot"], 
    # Scatter-specific parameters:
    markersize = 8, 
    color = :firebrick, 
    label = L"\alpha_{\mathrm{ov}}"
)
 ylims!(ax, 0, 4)
f
save("/home/ritavash/Desktop/Resources/convection_results/alpha_overshoot/16_history.png", f)
##
history = StellarModels.get_history_dataframe_from_hdf5("history.hdf5")
 history[!, "alpha_overshoot"][500]
#  history[!,"i_ov"][500]
##
using HDF5
using DataFrames# Assuming history is a Julia DataFrame

# Define the file name
filename = "history_32M.hdf5"

# Open the file in write mode ('w' will overwrite if it exists)
# Use the 'do' block to ensure the file is closed automatically
h5open(filename, "w") do file
    
    # Optional: Create a top-level group (e.g., "History") for organization
    group = create_group(file, "History")
    
    # Loop over all columns in the DataFrame
    for colname in names(history)
        # Write the column data to a dataset within the group
        group[colname] = history[!, colname]
    end
end

println("Successfully saved 'history' DataFrame to disk as $filename.")
##
### Perform some cleanup
#=
Internally we want to prevent storing any of the hdf5 files into our git repos, so I remove them. You can also take
advantage of `julia` as a scripting language to post-process your simulation output in a similar way.
=#
rm("history.hdf5")
rm("profiles.hdf5")
rm("example_options.toml")



