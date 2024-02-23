module Plotting

using GLMakie, LaTeXStrings, MathTeXEngine, Jems.StellarModels, Jems.DualSupport

const colors = Iterators.cycle([:red, :blue, :green])
const label_dict = Dict("mass" => L"m / M_\odot",
                        "zone" => L"\mathrm{zone}",
                        "dm" => L"dm / \mathrm{g}",

                        "dt" => L"dt / \mathrm{s}",
                        "star_age" => L"\mathrm{age} / \mathrm{year}",

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
                        "Y_center" => L"Y_\mathrm{center}"

                        )

include("Init.jl")
include("HRD.jl")
include("Profile.jl")
include("History.jl")

"""
    update_plots!(sm::StellarModel)

Updates all plots currently being displayed, by collecting appropriate data and notifying observables
"""
function update_plotting!(sm::StellarModel)
    if (sm.props.model_number % sm.opt.plotting.data_interval == 0)
        for plot in sm.plt.plots
            if plot.type == :HR
                update_HR_plot!(plot, sm.props)
            elseif plot.type == :profile
                update_profile_plot!(plot, sm)  # these cannot be loaded from props, bc they use the IO functions.
            elseif plot.type == :history
                update_history_plot!(plot, sm)
            end
        end
    end
    if (sm.props.model_number % sm.opt.plotting.plotting_interval == 0)
        for plot in sm.plt.plots
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
function end_of_evolution(sm::StellarModel)
    if sm.opt.plotting.wait_at_termination
        GLMakie.wait(sm.plt.scr)
    end
end

end  # end module Plotting
