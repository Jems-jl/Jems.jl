
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
# net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])


net = NuclearNetwork([:H1, :D2, :He3, :He4,
                      :Be7, :Li7, :B8
                      # :C12,   :C13,   
                      # :N13,   :N14,   :N15,
                      # :O14,   :O15,   :O16,   :O17,   :O18,   
                      # :F17,   :F18,   :F19,   
                      ],
                     [
                      # PP I
                    (:jina_rates, :H1_H1_to_D2_betplus_w_x_0),
                    (:jina_rates, :H1_H1_to_D2_xxec_w_x_0),
                    (:jina_rates, :H1_D2_to_He3_de04_n_x_0),
                    (:jina_rates, :H1_D2_to_He3_de04_x_x_0),
                    (:jina_rates, :He3_He3_to_H1_H1_He4_nacr_n_x_0),
                    # PP II
                    (:jina_rates, :He4_He3_to_Be7_cd08_n_x_0),
                    (:jina_rates, :He4_He3_to_Be7_cd08_n_x_1),
                    (:jina_rates, :Be7_to_Li7_xxec_w_x_0),
                    (:jina_rates, :H1_Li7_to_He4_He4_de04_x_x_0),
                    (:jina_rates, :H1_Li7_to_He4_He4_de04_r_x_0),
                    (:jina_rates, :H1_Li7_to_He4_He4_de04_x_x_1),
                    (:jina_rates, :H1_Li7_to_He4_He4_de04_r_x_1),
                    # PP III
                    (:jina_rates, :H1_Be7_to_B8_nacr_r_x_0),
                    (:jina_rates, :H1_Be7_to_B8_nacr_n_x_0),
                    (:jina_rates, :B8_to_He4_He4_wc12_w_x_0),
                    # PP IV
                    # (:jina_rates, :H1_He3_to_He4_betplus_w_x_0),

                     ])

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
          max_model_number = 1000
          max_center_T = 1e9

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
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0003, Chem.abundance_lists[:ASG_09], 
                                            1 * MSUN, 100 * RSUN; initial_dt=10 * SECYEAR)
@time Evolution.do_evolution_loop!(sm);


##

#=

Reached maximum model number
623.997675 seconds (153.16 M allocations: 10.026 GiB, 0.35% gc time, 2.77% compilation time: <1% of which was recompilation)


=#