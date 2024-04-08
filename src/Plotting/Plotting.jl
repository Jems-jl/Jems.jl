module Plotting

using GLMakie, LaTeXStrings, MathTeXEngine, Jems.StellarModels, Jems.DualSupport

const colors = Iterators.cycle(Makie.wong_colors())
const label_dict = Dict("mass" => L"m / M_\odot",
                        "zone" => L"\mathrm{zone}",
                        "dm" => L"dm / \mathrm{g}",
                        "model_number" => L"\mathrm{model number}",

                        "dt" => L"dt / \mathrm{s}",
                        "age" => L"\mathrm{age} / \mathrm{year}",

                        "log10_P" => L"\log_{10}(P / \mathrm{dyne})",
                        "log10_ρ" => L"\log_{10}(ρ / \mathrm{(g / cm^3)})",
                        "log10_T" => L"\log_{10}(T / \mathrm{K})",
                        "luminosity" => L"L / L_\odot",
                        "log10_r" => L"\log_{10}(r / \mathrm{cm})",
                        "X" => L"X",
                        "Y" => L"Y",

                        "P_surf" => L"P_\mathrm{surf} / \mathrm{dyne}",
                        "ρ_surf" => L"\rho_\mathrm{surf} / \mathrm{(g / cm^3})",
                        "T_surf" => L"T_\mathrm{surf} / \mathrm{K}",
                        "L_surf" => L"L_\mathrm{surf} / L_\odot",
                        "R_surf" => L"R_\mathrm{Surf} / R_\odot",
                        "X_surf" => L"X_\mathrm{surf}",
                        "Y_surf" => L"Y_\mathrm{surf}",

                        "P_center" => L"P_\mathrm{center} / \mathrm{dyne}",
                        "ρ_center" => L"\rho_\mathrm{center} / \mathrm{(g / cm^3})",
                        "T_center" => L"T_\mathrm{center} / \mathrm{K}",
                        "X_center" => L"X_\mathrm{center}",
                        "Y_center" => L"Y_\mathrm{center}",

                        "T" => L"T / \mathrm{K}",


                        "n" => L"n",
                        "H1" => L"^1H", "D2" => L"^2H",
                        "He3" => L"^3He", "He4" => L"^4He",
                        "Li7" => L"^7Li", "Be7" => L"^7Be",
                        "B8" => L"^8B"
                        )

include("Init.jl")
include("HRD.jl")
include("Profile.jl")
include("History.jl")

"""
    update_plots!(sm::StellarModel)

Updates all plots currently being displayed, by collecting appropriate data and notifying observables
"""
function update_plotting!(m::AbstractModel)
    # get new data
    if (m.props.model_number % m.opt.plotting.data_interval == 0)
        for plot in m.plt.plots
            if plot.type == :HR
                update_HR_plot!(plot, m.props)
            elseif plot.type == :profile
                update_profile_plot!(plot, m)  # these cannot be loaded from props, bc they use the IO functions.
            elseif plot.type == :history
                update_history_plot!(plot, m)
            end
        end
    end
    if (m.props.model_number % m.opt.plotting.plotting_interval == 0)
        for plot in m.plt.plots
            for xobs in values(plot.x_obs)
                notify(xobs)  # notifying only the x observables should replot everything
            end
            try
                autolimits!(plot.ax)
            catch  # catches Float32 errors in Makie
            end  # empty blocks allowed in julia!
            if !isnothing(plot.alt_ax)
                try
                    autolimits!(plot.alt_ax)
                catch
                end
            end
        end
    end
end

"""
    end_of_evolution(sm::StellarModel)

Perform end of evolution actions
"""
function end_of_evolution(m::AbstractModel)
    if m.opt.plotting.wait_at_termination
        GLMakie.wait(m.plt.scr)
    end
end

end  # end module Plotting
