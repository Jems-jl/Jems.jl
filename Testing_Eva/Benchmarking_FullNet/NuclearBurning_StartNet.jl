# NuclearBurning.jl

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

##

varnames = [:lnρ, :lnT, :lnr, :lum]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4], [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, structure_equations, nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);

##

open("example_options.toml", "w") do file
    write(file,
          """
          [remesh]
          do_remesh = true

          [solver]
          newton_max_iter_first_step = 1000
          initial_model_scale_max_correction = 0.5
          newton_max_iter = 30
          scale_max_correction = 0.1
          report_solver_progress = false

          [timestep]
          dt_max_increase = 1.5
          delta_R_limit = 0.01
          delta_Tc_limit = 0.01
          delta_Xc_limit = 0.005

          [termination]
          max_model_number = 2000
          max_center_T = 1e8

          [plotting]
          do_plotting = false
          wait_at_termination = false
          plotting_interval = 1

          window_specs = ["HR", "Kippenhahn", "profile", "TRhoProfile"]
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

          max_log_eps = 5.0

          [io]
          profile_interval = 50
          terminal_header_interval = 100
          terminal_info_interval = 100

          """)
end

StellarModels.set_options!(sm.opt, "./example_options.toml")
rm(sm.opt.io.hdf5_history_filename; force=true)
rm(sm.opt.io.hdf5_profile_filename; force=true)

n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], 
                                            1 * MSUN, 100 * RSUN; initial_dt=10 * SECYEAR)
@time Evolution.do_evolution_loop!(sm);

#=

TOY RATES

Reached maximum central temperature
149.255672 seconds (131.71 M allocations: 5.784 GiB, 1.83% gc time, 9.82% compilation time)

Reached maximum central temperature
151.596157 seconds (131.71 M allocations: 5.784 GiB, 1.90% gc time, 10.00% compilation time)

Reached maximum central temperature
159.586236 seconds (131.71 M allocations: 5.784 GiB, 1.96% gc time, 10.74% compilation time)

With plotting:

Reached maximum central temperature
464.764778 seconds (770.49 M allocations: 54.665 GiB, 2.66% gc time, 7.97% compilation time: 4% of which was recompilation)

=#

##

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
sm = StellarModel(varnames, structure_equations, nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);

##

open("example_options.toml", "w") do file
    write(file,
          """
          [remesh]
          do_remesh = true

          [solver]
          newton_max_iter_first_step = 1000
          initial_model_scale_max_correction = 0.5
          newton_max_iter = 30
          scale_max_correction = 0.1
          report_solver_progress = false

          [timestep]
          dt_max_increase = 1.5
          delta_R_limit = 0.01
          delta_Tc_limit = 0.01
          delta_Xc_limit = 0.005

          [termination]
          max_model_number = 2000
          max_center_T = 1e8

          [plotting]
          do_plotting = false
          wait_at_termination = false
          plotting_interval = 1

          window_specs = ["HR", "Kippenhahn", "profile", "TRhoProfile"]
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

          max_log_eps = 5.0

          [io]
          profile_interval = 50
          terminal_header_interval = 100
          terminal_info_interval = 100

          """)
end

StellarModels.set_options!(sm.opt, "./example_options.toml")
rm(sm.opt.io.hdf5_history_filename; force=true)
rm(sm.opt.io.hdf5_profile_filename; force=true)

n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], 
                                            1 * MSUN, 100 * RSUN; initial_dt=10 * SECYEAR)
@time Evolution.do_evolution_loop!(sm);


#=

Reached maximum central temperature
95.730732 seconds (86.59 M allocations: 3.263 GiB, 1.35% gc time, 18.24% compilation time)

Reached maximum central temperature
92.267297 seconds (86.72 M allocations: 3.272 GiB, 1.33% gc time, 18.68% compilation time)

Reached maximum central temperature
92.956283 seconds (86.75 M allocations: 3.272 GiB, 1.28% gc time, 19.99% compilation time)

Reached maximum central temperature
96.392549 seconds (86.72 M allocations: 3.272 GiB, 0.90% gc time, 16.86% compilation time)

 =#