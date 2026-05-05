using BenchmarkTools
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


net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
nz = 2000
nextra = 100
eos = EOS.IdealEOS(true)
op_es = Opacity.SimpleElectronScatteringOpacity() #for electron scattering 
"""
For using the opacity tables, you need to use the mesa opacity tables, with the path to directory and the identifier 
at the start. And then you need to use them with the transition region (in the same way as shown below )
"""
# low_T_collection = OpacityTableCollector("/Users/rdbnath/Documents/share_oplib_type1_tables/kap_low_alt/", "lowT_fa05_gs98") 
# high_T_collection = OpacityTableCollector("/Users/rdbnath/Documents/share_oplib_type1_tables/kap_data/","oplib_agss09" ) 
# op_oplib = CompositeOpacity(low_T_collection, high_T_collection, 3.8, 4.2)
turbulence = Turbulence.BasicMLT(2.0)
sm = StellarModel(TDCEquationSet(), nz, nextra, net, eos, op_oplib, turbulence);

##
open("example_options.toml", "w") do file
    write(file,
          """
          [remesh]
          do_remesh = true

          [solver]
          newton_max_iter_first_step = 1000
          initial_model_scale_max_correction = 0.2
          newton_max_iter = 200
          scale_max_correction = 0.1
          solver_progress_iter = 1
          relative_correction_tolerance = 1e12
          maximum_residual_tolerance = 1e-2
          use_preconditioning = true

          [timestep]
          dt_max_increase = 1.5
          delta_R_limit = 0.01
          delta_Tc_limit = 0.01
          delta_Xc_limit = 0.005

          [termination]
          max_model_number = 20000
          max_center_T = 1e8

          [io]
          profile_interval = 1
          terminal_header_interval = 100
          terminal_info_interval = 100
          profile_values = ["zone", "mass", "dm", "log10_rho", "log10_r", "log10_P", "log10_T", "luminosity",
                                      "X", "Y"]

          """)
end
StellarModels.set_options!(sm.opt, "./example_options.toml")
rm(sm.opt.io.hdf5_history_filename; force=true)
rm(sm.opt.io.hdf5_profile_filename; force=true)

##
#Configure live plots. To turn off one can use `plotter = Plotting.NullPlotter()`
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
#set initial condition and run model
n = 1.5
pms_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], 
                                            5 * MSUN, 100 * RSUN; initial_dt=0.01 * SECYEAR)
@time Evolution.do_evolution_loop!(sm, plotter=plotter);
 
##
