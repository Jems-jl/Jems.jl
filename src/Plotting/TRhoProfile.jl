"""
    function create_T_ρ_observables!(plot::StellarModels.JemsPlot, props::StellarModelProperties)

Creates relevant Tρ observables for this `plot` given the StellarModelProperties `props`
"""
function create_T_ρ_observables!(plot::StellarModels.JemsPlot, props::StellarModelProperties)
    rhos = get_value.(props.lnρ[1:(props.nz)]) .* log10(ℯ)
    ts = get_value.(props.lnT[1:props.nz]) .* log10(ℯ)
    plot.x_obs[:log_ρ] = Observable{Vector{Float64}}(rhos)
    plot.y_obs[:log_T] = Observable{Vector{Float64}}(ts)
    plot.other_obs[:colors] = 
            Observable{Vector{<:GLMakie.Colorant}}(mixing_colors[get.(Ref(mixing_map), props.mixing_type[1:(props.nz)],
                                                           missing)])
end


function _read_trho_data(file)
    ρs::Vector{Float64} = []
    Ts::Vector{Float64} = []
    open(file) do io
        for line in eachline(io)
            ρ, T = split(line)
            push!(ρs, parse(Float64, ρ))
            push!(Ts, parse(Float64, T))
        end
    end
    return ρs, Ts
end


"""
    function make_T_ρ_plot!(ax::Axis, ρ::Observable, T::Observable; line_kwargs=Dict())

Plots a line for the Tρ profile, along with burning lines, degeneracy line and Pgas ≈ Prad line. 
"""
function make_T_ρ_plot!(ax::Axis, ρ::Observable, T::Observable, mixing::Observable, line_kwargs=Dict())
    ax.xlabel = label_dict["log10_ρ"]
    ax.ylabel = label_dict["log10_T"]
    lines!(ax, ρ, T, line_kwargs..., color=mixing)
    h_burn_ρ, h_burn_T = _read_trho_data(pkgdir(Plotting, "data/PlotData", "hydrogen_burn.data"))
    lines!(ax, h_burn_ρ, h_burn_T, color=:gray, linestyle=:dash, linewidth=2)
    he_burn_ρ, he_burn_T = _read_trho_data(pkgdir(Plotting, "data/PlotData", "helium_burn.data"))
    lines!(ax, he_burn_ρ, he_burn_T, color=:gray, linestyle=:dash, linewidth=2)
    e_degen_ρ, e_degen_T = _read_trho_data(pkgdir(Plotting, "data/PlotData", "psi4.data"))
    lines!(ax, e_degen_ρ, e_degen_T, color=:gray, linestyle=:dash, linewidth=2)
    pgas_ρ = [-8, 5]
    pgas_T = log10(3.2e7) .+ (pgas_ρ .- log10(0.7e0))./3.0
    lines!(ax, pgas_ρ, pgas_T, color=:gray, linestyle=:dash, linewidth=2)
    nom = LineElement(color=mixing_colors[mixing_map[:no_mixing]])
    conv = LineElement(color=mixing_colors[mixing_map[:convection]])
    axislegend(ax, [nom, conv], ["no mixing", "convection"],position=:lt)
end

"""
    function update_T_ρ_plot!(plot::StellarModels.JemsPlot, props::StellarModelProperties)

updates the observables of this Tρ `plot` with relevant data of the stellar model properties `props`.
"""
function update_T_ρ_plot!(plot::StellarModels.JemsPlot, props::StellarModelProperties)
    plot.x_obs[:log_ρ].val = get_value.(props.lnρ[1:(props.nz)]) .* log10(ℯ)
    plot.y_obs[:log_T].val = get_value.(props.lnT[1:(props.nz)]) .* log10(ℯ)
    plot.other_obs[:colors].val = mixing_colors[get.(Ref(mixing_map), props.mixing_type[1:(props.nz)],
                                                           missing)]
end

