
using BenchmarkTools
using Jems.Chem
using Jems.Constants
using Jems.EOS
using Jems.Opacity
using Jems.NuclearNetworks
using Jems.StellarModels
using Jems.Evolution
using Jems.ReactionRates


##

varnames = [:lnρ, :lnT, :lnr, :lum]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1,:He4, :C12, :N14, :O16, :Ne20], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno), (:kipp_rates, :kipp_3alpha)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(false)
opacity = Opacity.SimpleElectronScatteringOpacity()
sm = StellarModel(varnames, structure_equations, nz, nextra,
                  remesh_split_functions, net, eos, opacity);

##

open("example_options.toml", "w") do file
    write(file,
          """
          [solver]
          newton_max_iter_first_step = 1000
          newton_max_iter = 200

          [timestep]
          dt_max_increase = 10.0
          delta_R_limit = 0.02
          delta_Tc_limit = 0.02

          [termination]
          max_model_number = 2000
          max_center_T = 1e9

          [io]
          profile_interval = 50
          """)
end
StellarModels.set_options!(sm.opt, "./example_options.toml")
rm(sm.opt.io.hdf5_history_filename; force=true)
rm(sm.opt.io.hdf5_profile_filename; force=true)
StellarModels.n1_polytrope_initial_condition!(sm, 10MSUN, 100 * RSUN; initial_dt=1000 * SECYEAR)

@time Evolution.do_evolution_loop(sm);

# max_center_T = 1e9 --> this T is the He burning temperature

##

using CairoMakie, LaTeXStrings, MathTeXEngine
basic_theme = Theme(
                    fonts = (regular = texfont(:text), bold = texfont(:bold),
                    italic = texfont(:italic), bold_italic = texfont(:bolditalic)),
                    fontsize=30, resolution=(1000, 750), linewidth=7,
                    Axis=(xlabelsize=40, ylabelsize=40, titlesize=40, xgridvisible=false, ygridvisible=false,
                          spinewidth=2.5, xminorticksvisible=true, yminorticksvisible=true, xtickalign=1, ytickalign=1,
                          xminortickalign=1, yminortickalign=1, xticksize=14, xtickwidth=2.5, yticksize=14,
                          ytickwidth=2.5, xminorticksize=7, xminortickwidth=2.5, yminorticksize=7, yminortickwidth=2.5,
                          xticklabelsize=35, yticklabelsize=35, xticksmirrored=true, yticksmirrored=true),
                    Legend=(patchsize=(70, 10), framevisible=false, patchlabelgap=20, rowgap=10))
set_theme!(basic_theme)

##

f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(T_\mathrm{eff}/[K])", ylabel=L"\log_{10}(L/L_\odot)", xreversed=true)
history = StellarModels.get_history_dataframe_from_hdf5("history.hdf5")
scatter!(ax, log10.(history[!, "T_surf"]), log10.(history[!, "L_surf"]))
f
# save("HR_1.4M_140R.png", f, px_per_unit = 2)

##

f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(T_\mathrm{eff}/[K])", ylabel=L"\log_{10}(L/L_\odot)", xreversed=true)
history = StellarModels.get_history_dataframe_from_hdf5("history.hdf5")
scatter!(ax, (history[!, "X_center"]), 1 .- history[!, "X_center"] .- history[!, "Y_center"])
f

##

f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(T_\mathrm{eff}/[K])", ylabel=L"\log_{10}(L/L_\odot)", xreversed=true)
history = StellarModels.get_history_dataframe_from_hdf5("history.hdf5")
scatter!(ax, (history[!, "model_number"]), history[!, "T_center"])
f

##

He4 = StellarModels.history_get_ind_vars_edge_value(sm, :He4, :center)
H1 = StellarModels.history_get_ind_vars_edge_value(sm, :H1, :center)
1 - H1 - He4

##

f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(T_\mathrm{eff}/[K])", ylabel=L"\log_{10}(L/L_\odot)")

profile_names = StellarModels.get_profile_names_from_hdf5("profiles.hdf5")

profile = StellarModels.get_profile_dataframe_from_hdf5("profiles.hdf5", "0000000600")
mass = profile[!, "mass"]
X = profile[!, "X"]

lines!(ax, mass, X; label="real profile")
f
 