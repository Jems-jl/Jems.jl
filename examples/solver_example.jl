#=
# solver_example.jl

This example demonstrates the evolution solver interface:

  * configuring `NewtonSolver` through `DefaultEvolver`;
  * retrying failed steps with smaller timesteps;
  * falling back to the original mesh after remeshed attempts fail;
  * keeping remeshing disabled for a configurable number of accepted steps;
  * adding a user-defined termination condition;
  * overriding the default timestep controller; and
  * running hooks before and after stellar-model property evaluations.
=#

using Jems.Chem
using Jems.Constants
using Jems.EOS
using Jems.Opacity
using Jems.NuclearNetworks
using Jems.Turbulence
using Jems.StellarModels
using Jems.Evolution
using Jems.Plotting
using Jems.Interpolations

##
# Download and extract the MESA EOS and opacity data when they are not already
# available. This is the same data preparation used by opacity_example.jl.
if !isdir("MESA_data")
    mkdir("MESA_data")
    download("https://zenodo.org/records/19722306/files/mesa-26.04.1.zip",
             "MESA_data/mesa-26.04.1.zip")

    cd("MESA_data")
    run(`unzip -q mesa-26.04.1.zip`)

    mv("mesa-26.04.1/kap/kap_data.tar.xz", "./kap_data.tar.xz")
    run(`tar -xJf kap_data.tar.xz`)
    rm("kap_data.tar.xz")

    mkdir("eosFreeEOS_data_filtered")
    mv("mesa-26.04.1/eos/eosFreeEOS_data.tar.xz", "./eos_data.tar.xz")
    run(`tar -xJf eos_data.tar.xz`)
    rm("eos_data.tar.xz")

    extracted_eos_dir = "eosFreeEOS_data"
    target_pattern = r"^mesa-FreeEOS_\d+z\d+x\.data$"
    for file in readdir(extracted_eos_dir)
        if occursin(target_pattern, file)
            mv(joinpath(extracted_eos_dir, file),
               joinpath("eosFreeEOS_data_filtered", file))
        end
    end

    rm("mesa-26.04.1", recursive=true)
    rm("eosFreeEOS_data", recursive=true)
    mv("eosFreeEOS_data_filtered", "eosFreeEOS_data")
    rm("mesa-26.04.1.zip")
    cd("..")
end

##
# Microphysics and stellar model.
net = merge_nuclear_networks([
    NuclearNetworks.networks[:basic_He_burn],
    NuclearNetworks.networks[:basic_CO_burn],
    NuclearNetworks.networks[:JINA_PPI],
    NuclearNetworks.networks[:JINA_PPII],
    NuclearNetworks.networks[:JINA_PPIII],
    NuclearNetworks.networks[:JINA_CNOI],
    NuclearNetworks.networks[:JINA_CNOII],
])

nz = 1000
nextra = 100
eos_table = EOSTableCollector( "MESA_data/eosFreeEOS_data"; include_radiation=true,interp_type=:bilinear,)
low_T_collection = OpacityTableCollector( "MESA_data/kap_data","lowT_fa05_gs98",)
high_T_collection = OpacityTableCollector( "MESA_data/kap_data","oplib_agss09",)
opacity = CompositeOpacity(low_T_collection, high_T_collection, 3.8, 4.2)
turbulence = Turbulence.BasicMLT(2.0)

sm = StellarModel( StellarModels.DefaultStellarEquationSet(),nz, nextra, net, eos_table, opacity, turbulence,);

##
# Standard model options remain in the TOML file. Recovery behavior that is
# part of the evolution algorithm is configured on DefaultEvolver below.
open("solver_example_options.toml", "w") do file
    write(file,
          """
          [remesh]
          do_remesh = true

          [solver]
          newton_max_iter_first_step = 1000
          initial_model_scale_max_correction = 0.2
          newton_max_iter = 100
          scale_max_correction = 0.1
          report_solver_progress = true
          solver_progress_iter = 10
          relative_correction_tolerance = 1e12
          maximum_residual_tolerance = 1e-2
          use_preconditioning = true

          [physics]
          freeze_D_turb_during_step = true

          [timestep]
          dt_max_increase = 1.5
          dt_max_decrease = 0.5
          dt_retry_decrease = 0.5
          delta_R_limit = 0.01
          delta_Tc_limit = 0.01
          delta_Xc_limit = 0.005

          [termination]
          # This is deliberately larger than HOOK_MODEL_LIMIT below so that
          # the custom termination hook is what stops this example.
          max_model_number = 1000
          max_center_T = 1e40

          [io]
          hdf5_history_filename = "solver_history.hdf5"
          hdf5_profile_filename = "solver_profiles.hdf5"
          profile_interval = 25
          terminal_header_interval = 25
          terminal_info_interval = 25
          profile_values = ["zone", "mass", "dm", "log10_rho", "log10_r","log10_P", "log10_T", "luminosity", "X", "Y"]
          """)
end

StellarModels.set_options!(sm.opt, "./solver_example_options.toml")
rm(sm.opt.io.hdf5_history_filename; force=true)
rm(sm.opt.io.hdf5_profile_filename; force=true)


##
# Hook configuration example for timestep control, termination control 
const HOOK_MODEL_LIMIT = 200
const HOOK_TIMESTEP_CAP = 1e6 * SECYEAR

"""Stop this example through the user termination hook."""
function Evolution.hook_termination_condition( sm::StellarModel, ::Evolution.DefaultEvolver,)
    sm.props.model_number >= HOOK_MODEL_LIMIT && return :example_model_limit
    return nothing
end

"""Cap the default timestep after it has been computed by Jems."""
function Evolution.hook_get_dt_next(sm::StellarModel, ::Evolution.DefaultEvolver,)
    sm.props.dt_next > HOOK_TIMESTEP_CAP && return HOOK_TIMESTEP_CAP
    return nothing
end


##
# The evolver selects the nonlinear solver and controls recovery from failed
# steps. If remeshed attempts are exhausted, the current step is retried on the
# original mesh. After it succeeds, the following five accepted steps also use
# the original mesh before remeshing is enabled again.
evolver = Evolution.DefaultEvolver( solver=Evolution.NewtonSolver(), max_timestep_retries=10, skip_remesh_steps=5,)

##
# Configure live plots. Replace this block with
# `plotter = Plotting.NullPlotter()` to run without a graphical window.
using GLMakie
GLMakie.activate!()
set_theme!(Plotting.basic_theme())
f = Figure(size=(1400,750))
plots = [Plotting.HRPlot(f[1,1]),
         Plotting.TRhoProfile(f[1,2]),
         Plotting.KippenLine(f[2,1], xaxis=:time, time_units=:Gyr),
         Plotting.AbundancePlot(f[2,2],net,log_yscale=true, ymin=1e-3),
         Plotting.HistoryPlot(f[1,3], sm, x_name="age", y_name="X_center", othery_name="Y_center", link_yaxes=true),
         Plotting.ProfilePlot(f[2,3], sm, x_name="mass", y_name="log10_rho", othery_name="log10_T")]
plotter = Plotting.Plotter(fig=f,plots=plots)


##
# Initial model and evolution.
n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], 
                                             1 * MSUN, 20 * RSUN; initial_dt=0.01 * SECYEAR)

@time Evolution.do_evolution_loop!(sm, plotter=plotter);
