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
using Jems.DualSupport
using ForwardDiff
tag_external = ForwardDiff.Tag{:external, nothing}
nbmodmax = 4000
ForwardDiff.tagcount(tag_external); #this function is necessary to order the tags
open("example_options.toml", "w") do file
    write(file,
          """
          [remesh]
          do_remesh = true

          [solver]
          newton_max_iter_first_step = 1000
          newton_max_iter = 30
          initial_model_scale_max_correction = 0.1
          scale_max_correction = 0.1
          report_solver_progress = false

          [timestep]
          dt_max_increase = 1.5
          delta_R_limit = 0.01
          delta_Tc_limit = 0.01
          delta_Xc_limit = 0.002

          [termination]
          max_model_number = $nbmodmax
          max_center_T = 1e8

          [plotting]
          do_plotting = false
          wait_at_termination = false
          plotting_interval = 1

          window_specs = ["HR", "TRho", "profile", "history"]
          window_layouts = [[1, 1],  # arrangement of plots
                            [1, 2],
                            [2, 1],
                            [2, 2]
                            ]

          profile_xaxis = 'mass'
          profile_yaxes = ['log10_T']
          profile_alt_yaxes = ['X','Y']

          history_xaxis = 'star_age'
          history_yaxes = ['R_surf']
          history_alt_yaxes = ['T_center']

          [io]
          history_interval = 1
          profile_interval = 50
          terminal_header_interval = 50
          terminal_info_interval = 10

          """)
end
##
logM_range = [0.15,0.25]
X = 0.7381
overwrite = false
## Model creation, as usual
for logM in logM_range
    M = 10^logM
    println("############################################################################################")
    println("STARTING NEW logM = $logM , M = $M ###########################################")
    history_path = "Jems.jl/DualRuns/DualGrid/" * "logM_" * string(logM) * "_" * "X_" * string(X) * "_" * ".history.hdf5"
    profile_path = "Jems.jl/DualRuns/DualGrid/" * "logM_" * string(logM) * "_" * "X_" * string(X) * "_" * ".profiles.hdf5"
    #history_path = "DualRuns/DualGrid/" * "logM_" * string(logM) * "_" * "X_" * string(X) * "_" * ".history.hdf5"
    #profile_path = "DualRuns/DualGrid/" * "logM_" * string(logM) * "_" * "X_" * string(X) * "_" * ".profiles.hdf5"
    
    @show history_path
    @show profile_path
    if isfile(history_path)
        if overwrite == true
            println("History file for mass = $M / logM = $logM already exists: OVERWRITING")
        else
            println("History file for mass = $M / logM = $logM already exists: SKIPPING")
            continue
        end
    end


    println("Create StellarModel ############################")
    varnames = [:lnρ, :lnT, :lnr, :lum]
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
    number_of_partials = 5 #number of partials we want to keep track of
    dual_type = ForwardDiff.Dual{tag_external,Float64,number_of_partials}
    sm = StellarModel(varnames, structure_equations, nz, nextra, remesh_split_functions, net, eos, opacity, turbulence, number_type = dual_type);
    ##
    println("Initialize StellarModel ############################")
    n = 3
    #define dual input numbers, all partial derivatives are with respect to the mass, give a '1.0' to indicate the order of the partials
    logM_dual      = ForwardDiff.Dual{tag_external}(logM,     1.0,0.0,0.0,0.0,0.0)
    mass_dual      = MSUN*10^logM_dual
    X_dual         = ForwardDiff.Dual{tag_external}(X,  0.0,1.0,0.0,0.0,0.0)
    Z_dual         = ForwardDiff.Dual{tag_external}(0.0134,  0.0,0.0,1.0,0.0,0.0)
    #Dfraction_dual = ForwardDiff.Dual{tag_external}(0.0,     0.0,0.0,0.0,1.0,0.0)
    Dfraction_dual = ForwardDiff.Dual{tag_external}(0.000312,0.0,0.0,0.0,1.0,0.0)
    R_dual         = ForwardDiff.Dual{tag_external}(100*RSUN,0.0,0.0,0.0,0.0,1.0)
    StellarModels.n_polytrope_initial_condition!(n, sm, nz, X_dual,Z_dual,Dfraction_dual,Chem.abundance_lists[:ASG_09],mass_dual, R_dual; initial_dt=10 * SECYEAR)
    StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
    Evolution.cycle_props!(sm);
    StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
    StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
    StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
    StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
    StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)
    print("Initialize & Evolve StellarModel ############################################################")
    StellarModels.set_options!(sm.opt, "./example_options.toml")
    rm(sm.opt.io.hdf5_history_filename; force=true)
    rm(sm.opt.io.hdf5_profile_filename; force=true)
    StellarModels.n_polytrope_initial_condition!(n, sm, nz, X_dual,Z_dual,Dfraction_dual,Chem.abundance_lists[:ASG_09],
                                                mass_dual, R_dual; initial_dt=10 * SECYEAR)
    @time Evolution.do_evolution_loop!(sm);
    cp("history.hdf5",  history_path, force = true)
    cp("profiles.hdf5", profile_path, force = true)
    println("RUN ENDED, SAVED HISTORY AND PROFILES")
    @show history_path
    @show profile_path
end
##
history_path = "Jems.jl/DualRuns/DualGrid/logM_0.03_X_0.7381_.history.hdf5"
profile_path = "Jems.jl/DualRuns/DualGrid/logM_0.03_X_0.7381_.profiles.hdf5"
cp("history.hdf5",  history_path, force = true)
cp("profiles.hdf5", profile_path, force = true)
println("RUN ENDED, SAVED HISTORY AND PROFILES")
@show history_path
@show profile_path
##
#=
###Accessing the Dual profiles

For the history: it can still be accessed as before, but now we also have access to the partials.
For the profiles: it can still be accessed as before, but we should take care not to bump into errors. 
Therefore i give a `value_names`, the string names of the original profiles, and `dual_names`, 
containing lists of the string names of the corresponding partial profiles.
=#
##  
using DataFrames
using HDF5
"""
get_dual_profile_dataframe_from_hdf5(hdf5_filename, value_name, partials_names)

Returns a DataFrame object built filled with Dual numbers.
"""
function get_partial_profile_dataframe_from_hdf5(hdf5_filename, value_name, partials_names)
    value_dataframe = StellarModels.get_profile_dataframe_from_hdf5(hdf5_filename, value_name)
    partial_dataframes = [StellarModels.get_profile_dataframe_from_hdf5(hdf5_filename, partial_name) for partial_name in partials_names]
    df_partial = ForwardDiff.Dual.(value_dataframe, partial_dataframes...)
    return df_partial
end
 #doing some bookkeeping stuff
profile_names = StellarModels.get_profile_names_from_hdf5("profiles.hdf5")#all profiles, regular profiles and dual profiles
value_names = [name for name in profile_names if !occursin("partial", name)]
partial_names_unpacked = [name for name in profile_names if occursin("partial", name)]
partial_names = [[partial_name for partial_name in partial_names_unpacked[lo:lo+number_of_partials-1] ] for lo in 1:number_of_partials:(length(partial_names_unpacked))] 


################################################################################## PROFILE OUTPUT
value_names #this list contains the profile names as before, i.e. just the values, nothing special
partial_names #this list contains lists with the corresponding partial names
i = 2
StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", value_names[i]) #access the ith profile with actual values, as before
bla = get_partial_profile_dataframe_from_hdf5("profiles.hdf5", value_names[i], partial_names[i]) #acces the ith profile, but now with Dual numbers, i.e. containg both the values and the partials  
#################################################################################

function get_dual_history_dataframe_from_hdf5(hdf5_filename)
    #This function used two functions that were originally defined for profile handling, but they come in handy here
    names = StellarModels.get_profile_names_from_hdf5(hdf5_filename)
    @show names
    history_value_name = names[1]
    history_partial_names = names[2:end]
    history_value = StellarModels.get_history_dataframe_from_hdf5(hdf5_filename)
    history_partials = [StellarModels.get_profile_dataframe_from_hdf5(hdf5_filename, name) for name in history_partial_names]
    return ForwardDiff.Dual.(history_value, history_partials...)
end

################################################################################# HISTORY OUTPUT
history = StellarModels.get_history_dataframe_from_hdf5("history.hdf5") #as before
get_dual_history_dataframe_from_hdf5("history.hdf5")#dataframe with history in dual numbers

#################################################################################





























##
#=
### Plotting with Makie

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

Copying, for clarity, the same example from NuclearBurning.jl
=#

f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(\rho/\mathrm{[g\;cm^{-3}]})", ylabel=L"\log_{10}(P/\mathrm{[dyn]})")

pname = Observable(value_names[end])#NOTE: we use value_names here!

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

f

##
record(f, "rho_P_evolution.gif", value_names[1:end]; framerate=2) do profile_name
    profile = StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", profile_name)
    log10_P.val = profile[!, "log10_P"]
    pname[] = profile_name
end

##

#=
### Check nuclear burning


=#

f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\mathrm{Mass}\;[M_\odot]", ylabel=L"X")

pname = Observable(value_names[end]) #NOTE that we used value_names here!

profile = @lift(StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", $pname))
mass = @lift($profile[!, "mass"])
X = Observable(rand(length(mass.val)))
model_number_str = @lift("model number=$(parse(Int,$pname))")

profile_line = lines!(ax, mass, X; label="real profile")
profile_text = text!(ax, 0.7, 0.0; text=model_number_str)
f
record(f, "X_evolution.gif", value_names[1:end]; framerate=2) do profile_name
    profile = StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", profile_name)
    X.val = profile[!, "X"]
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
history = StellarModels.get_history_dataframe_from_hdf5("history.hdf5") #NOTE that history file just works as before
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
