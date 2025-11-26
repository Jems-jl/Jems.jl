module Plotting

using Makie, LaTeXStrings, MathTeXEngine, Jems.DualSupport, Jems.Constants

include("Plotter.jl")
include("HRPlot.jl")
include("TRhoProfile.jl")

# using GLMakie, LaTeXStrings, MathTeXEngine, Jems.StellarModels, Jems.DualSupport, Jems.Constants
# 
# const colors = Iterators.cycle(Makie.wong_colors())
# const mixing_map = Dict(:no_mixing => 1,
#                         :convection => 2)
# mixing_colors = [RGBAf(0.5, 0.5, 0.5), RGBAf(0, 0, 1)]  # mixing colors in TRhoProfile diagrams
# kipp_mixing_colors = copy(mixing_colors)
# kipp_mixing_colors[1] = RGBAf(1, 1, 1)  # make no_mixing white in KippenLine diagram to avoid clutter
# burning_colors = cgrad(:linear_wyor_100_45_c55_n256)
# function burning_map(log_eps_nuc; min_log_eps=0.0, max_log_eps=15.0)  # map log eps nuc to interval [0.0, 1.0]
#     if log_eps_nuc < min_log_eps
#         return 0.0
#     elseif log_eps_nuc > max_log_eps
#         return 1.0
#     else
#         return log_eps_nuc / (max_log_eps - min_log_eps)
#     end
# end
# 
# const label_dict = Dict("mass" => L"m / M_\odot",
#                         "zone" => L"\mathrm{zone}",
#                         "dm" => L"dm / \mathrm{g}",
#                         "model_number" => "model number",
# 
#                         "dt" => L"dt / \mathrm{s}",
#                         "age" => L"\mathrm{age} / \mathrm{year}",
# 
#                         "log10_P" => L"\log_{10}(P / \mathrm{dyne})",
#                         "log10_ρ" => L"\log_{10}(ρ / \mathrm{(g / cm^3)})",
#                         "log10_T" => L"\log_{10}(T / \mathrm{K})",
#                         "luminosity" => L"L / L_\odot",
#                         "log10_r" => L"\log_{10}(r / \mathrm{cm})",
#                         "X" => L"X",
#                         "Y" => L"Y",
# 
#                         "P_surf" => L"P_\mathrm{surf} / \mathrm{dyne}",
#                         "ρ_surf" => L"\rho_\mathrm{surf} / \mathrm{(g / cm^3})",
#                         "T_surf" => L"T_\mathrm{surf} / \mathrm{K}",
#                         "L_surf" => L"L_\mathrm{surf} / L_\odot",
#                         "R_surf" => L"R_\mathrm{Surf} / R_\odot",
#                         "X_surf" => L"X_\mathrm{surf}",
#                         "Y_surf" => L"Y_\mathrm{surf}",
# 
#                         "P_center" => L"P_\mathrm{center} / \mathrm{dyne}",
#                         "ρ_center" => L"\rho_\mathrm{center} / \mathrm{(g / cm^3})",
#                         "T_center" => L"T_\mathrm{center} / \mathrm{K}",
#                         "X_center" => L"X_\mathrm{center}",
#                         "Y_center" => L"Y_\mathrm{center}",
# 
#                         "H1" => L"^1H", "D2" => L"^2H",
#                         "He3" => L"^3He", "He4" => L"^4He",
#                         "Li7" => L"^7Li", "Be7" => L"^7Be",
#                         "B8" => L"^8B"
#                         )
# 
# include("Init.jl")
# include("HRD.jl")
# include("Profile.jl")
# include("History.jl")
# include("TRhoProfile.jl")
# include("KippenLine.jl")
# 
# """
#     update_plots!(sm::StellarModel)
# 
# Updates all plots currently being displayed, by collecting appropriate data and notifying observables
# """
# function update_plotting!(m::AbstractModel)
#     # get new data
#     if (m.props.model_number % m.opt.plotting.data_interval == 0)
#         for plot in m.plt.plots
#             if plot.type == :HR
#                 update_HR_plot!(plot, m.props)
#             elseif plot.type == :TRhoProfile
#                 update_T_ρ_plot!(plot, m.props)
#             elseif plot.type == :profile
#                 update_profile_plot!(plot, m)  # these cannot be loaded from props, bc they use the IO functions.
#             elseif plot.type == :history
#                 update_history_plot!(plot, m)
#             elseif plot.type == :Kippenhahn
#                 update_Kipp_plot!(plot, m.props, m.opt.plotting)
#             end
#         end
#     end
#     if (m.props.model_number % m.opt.plotting.plotting_interval == 0)
#         for plot in m.plt.plots
#             for obs in values(plot.x_obs)
#                 notify(obs)  # notifying only the x observables should replot everything
#             end
#             if !isnothing(plot.other_obs)
#                 for obs in values(plot.other_obs)
#                     notify(obs)
#                 end
#             end
#             try
#                 autolimits!(plot.ax)
#             catch  # catches Float32 errors in Makie
#             end  # empty blocks allowed in julia!
#             if !isnothing(plot.alt_ax)
#                 try
#                     autolimits!(plot.alt_ax)
#                 catch
#                 end
#             end
#         end
#     end
# end
# 
# """
#     end_of_evolution(sm::StellarModel)
# 
# Perform end of evolution actions
# """
# function end_of_evolution(m::AbstractModel)
#     if m.opt.plotting.wait_at_termination
#         GLMakie.wait(m.plt.scr)
#     end
# end

end  # end module Plotting
