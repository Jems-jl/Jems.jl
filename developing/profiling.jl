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
using Profile
using PProf
using Jems.DualSupport
using ForwardDiff
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
n=3
StellarModels.n_polytrope_initial_condition!(n, sm, MSUN, 100 * RSUN; initial_dt=10 * SECYEAR)
Evolution.set_step_info!(sm, sm.esi)
Evolution.cycle_step_info!(sm);
Evolution.set_step_info!(sm, sm.ssi)

### Benchmarking
##
@benchmark begin
    my_update_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end


##
#=
To get an idea of how much a complete iteration of the solver takes, we need to benchmark
both the calculation of the Jacobian and the matrix solver. This is because the matrix solver
is destructive, as it uses the allocated Jacobian to store intermediate results. The time it takes
to run only the matrix solver can be determined by substracting the previous benchmark from this one.
=#

@benchmark begin
    StellarModels.update_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

##
Profile.clear()
Profile.@profile begin
    StellarModels.update_stellar_model_properties!(sm, sm.props)
    Evolution.eval_jacobian_eqs!(sm)
    Evolution.thomas_algorithm!(sm)
end
PProf.pprof()

##
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
          max_model_number = 2000
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
StellarModels.n_polytrope_initial_condition!(n, sm, 1 * MSUN, 100 * RSUN; initial_dt=1000 * SECYEAR)
@profview sm = Evolution.do_evolution_loop!(sm);

##
@benchmark StellarModels.write_data($sm)
