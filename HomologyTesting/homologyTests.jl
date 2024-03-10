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


#taking a new equation for the temperature, only including radiation
function equationT(sm::StellarModel, k::Int)
    lnT₀ = get_00_dual(sm.props.eos_res[k].lnT)
    if k == sm.props.nz  # atmosphere boundary condition
        L₀ = get_00_dual(sm.props.L[k]) * LSUN
        r₀ = exp(get_00_dual(sm.props.lnr[k]))
        return lnT₀ - log(L₀ / (BOLTZ_SIGMA * 4π * r₀^2)) / 4  # Eddington gray, ignoring radiation pressure term
    end
    r₀ = exp(get_00_dual(sm.props.lnr[k]))
    lnT₀ = get_00_dual(sm.props.lnT[k])
    lnT₊ = get_p1_dual(sm.props.lnT[k+1])

    Pface = exp(get_00_dual(sm.props.lnP_face[k]))
    Tface = exp(get_00_dual(sm.props.lnT_face[k]))

    ∇ᵣ = get_00_dual(sm.props.turb_res[k].∇ᵣ)
    return (Tface * (lnT₊ - lnT₀) / sm.props.dm[k] +
            CGRAV * sm.props.m[k] * Tface / (4π * r₀^4 * Pface) * ∇ᵣ) /
           (CGRAV * sm.props.m[k] * Tface / (4π * r₀^4 * Pface))

end

open("example_options.toml", "w") do file
    write(file,
          """
          [remesh]
          do_remesh = true

          [solver]
          newton_max_iter_first_step = 1000
          initial_model_scale_max_correction = 0.1
          newton_max_iter = 30
          scale_max_correction = 0.1
          report_solver_progress = false

          [timestep]
          dt_max_increase = 1.5
          delta_R_limit = 0.01
          delta_Tc_limit = 0.01
          delta_Xc_limit = 0.002

          [termination]
          max_model_number = 1400
          max_center_T = 1e8

          [plotting]
          do_plotting = true
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
          terminal_header_interval = 500
          terminal_info_interval = 200

          """)
end


##

logmassrange = (-1:0.1:1)
Xrange = (0.6:0.05:0.9)
#Xrange = (0.6:0.1:0.9)

luminosities = zeros(length(logmassrange), length(Xrange))
Z = 0.0134
Dfraction = 0.0
println(" BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN BEGIN ********************************")
for (i, logmass) in enumerate(logmassrange)
    println(" ---------------------------------------------------------------------------------")
    println(" ---------------------------------------------------------------------------------")
    mass = 10^logmass
    println("Mass = $mass, logM = $logmass")

    for (j,X) in enumerate(Xrange)

        #if history file already exists, skip
        path = "HomologyTesting/Histories/" * "toyRates_" * "logM"*string(logmass)*"_X"*string(X)*"_.hdf5"
        if isfile(path)
            println("history file for mass = $mass / logM = $logmass and X = $X already exists, skipping")
            continue
        end

        println(" ---------------------------------------------------------------------------------")
        println("Mass = $mass, logM = $logmass, X = $X")

        ## Model creation
        varnames = [:lnρ, :lnT, :lnr, :lum]
        structure_equations = [Evolution.equationHSE, Evolution.equationT,
                               Evolution.equationContinuity, Evolution.equationLuminosity]
        remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                                  StellarModels.split_lnT, StellarModels.split_xa]
        net = NuclearNetwork([:H1,:He4,:C12,:N14, :O16], [(:toy_rates, :toy_pp)])
        nz = 1000
        nextra = 1000
        eos = EOS.IdealEOS(false)
        opacity = Opacity.SimpleElectronScatteringOpacity()
        turbulence = Turbulence.BasicMLT(1.0)
        sm = StellarModel(varnames, structure_equations, nz, nextra, remesh_split_functions, net, eos, opacity, turbulence)
        StellarModels.set_options!(sm.opt, "./example_options.toml")
        rm(sm.opt.io.hdf5_history_filename; force=true)
        rm(sm.opt.io.hdf5_profile_filename; force=true)
        n = 3
        StellarModels.n_polytrope_initial_condition!(n,sm,nz,X,Z,Dfraction,Chem.abundance_lists[:ASG_09],mass*MSUN, mass^(3/2)*100 * RSUN;initial_dt=10 * SECYEAR)
        @time Evolution.do_evolution_loop!(sm);
        history = StellarModels.get_history_dataframe_from_hdf5("history.hdf5")
        
        cp("history.hdf5", path, force = true)
        println("saved history to $path")
        #luminosities[i,j] = history[!, "L_surf"][end]
        sleep(100)
    end
end
##
history
##
#save luminosities to .txt
using DelimitedFiles
writedlm("examples/luminosities.txt", luminosities)



##
5^(3/2) 
##


number = 5
number2 = 2.0
mystring = "Hello, world!"
#now add number to mystring
string(number)*mystring 
