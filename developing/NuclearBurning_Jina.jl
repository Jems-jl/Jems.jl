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
using Jems.Turbulence
using Jems.StellarModels
using Jems.Evolution
using Jems.ReactionRates

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
# net = NuclearNetwork([:H1,:He4,:C12,:N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])

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
                      # (:jina_rates, :H1_D2_to_He3_de04_x_x_0),
                      (:jina_rates, :He3_He3_to_H1_H1_He4_nacr_n_x_0),
                      # PP II
                      (:jina_rates, :He4_He3_to_Be7_cd08_n_x_0),
                      (:jina_rates, :He4_He3_to_Be7_cd08_n_x_1),
                      (:jina_rates, :Be7_to_Li7_xxec_w_x_0),
                      (:jina_rates, :H1_Li7_to_He4_He4_de04_x_x_0),
                      (:jina_rates, :H1_Li7_to_He4_He4_de04_r_x_0),
                      # (:jina_rates, :H1_Li7_to_He4_He4_de04_x_x_1),
                      # (:jina_rates, :H1_Li7_to_He4_He4_de04_r_x_1),
                      # PP III
                      (:jina_rates, :H1_Be7_to_B8_nacr_r_x_0),
                      (:jina_rates, :H1_Be7_to_B8_nacr_n_x_0),
                      (:jina_rates, :B8_to_He4_He4_wc12_w_x_0),
                      # PP IV
                      (:jina_rates, :H1_He3_to_He4_betplus_w_x_0)])

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
          do_remesh = false

          [solver]
          newton_max_iter_first_step = 1000
          initial_model_scale_max_correction = 0.2
          newton_max_iter = 50
          scale_max_correction = 0.2
          report_solver_progress = true
          solver_progress_iter = 50
          relative_correction_tolerance = 1e6
          maximum_residual_tolerance = 1e-6

          [timestep]
          dt_max_increase = 1.1
          delta_R_limit = 0.01
          delta_Tc_limit = 0.01
          delta_Xc_limit = 0.005

          [termination]
          max_model_number = 2000
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
          terminal_header_interval = 100
          terminal_info_interval = 100

          """)
end
StellarModels.set_options!(sm.opt, "./example_options.toml")
rm(sm.opt.io.hdf5_history_filename; force=true)
rm(sm.opt.io.hdf5_profile_filename; force=true)
n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0003, Chem.abundance_lists[:ASG_09],
                                             0.5 * MSUN, 50 * RSUN; initial_dt=1 * SECYEAR)

##
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
rm("history.hdf5")
rm("profiles.hdf5")
rm("example_options.toml")
