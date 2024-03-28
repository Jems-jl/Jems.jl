# """
#     function create_Kipp_observables!(plot::StellarModels.JemsPlot, props::StellarModelProperties)

# Creates relevant observables for this `plot` given the StellarModelProperties `props`
# """
# function create_Kipp_observables!(plot::StellarModels.JemsPlot, props::StellarModelProperties)
#     plot.y_obs[:mass] = Observable(Float64[])
#     push!(plot.y_obs[:mass][], props.mstar)
#     plot.other_obs[:mix_colors] = Observable(Vector{<:GLMakie.Colorant}[])
#     push!(plot.other_obs[:mix_colors][],
#           mixing_colors[get.(Ref(mixing_map), props.mixing_type[1:(props.nz)], missing)])
#     plot.other_obs[:burn] = Observable(Vector{<:GLMakie.Colorant}[])
#     push!(plot.other_obs[:burn][], props.eps_nuc[1:(props.nz)])
# end

"""
    function make_Kipp_plot!(ax::Axis)

Plots two lines per model in a Kippenhahn-like diagram (I dub it the KippenLine diagram).
One line contains info on the mixing state, the other on the burning regions.
"""
function make_Kipp_plot!(ax::Axis, model_number, mass, mixing, burn; line_kwargs=Dict())
    ax.xlabel = label_dict["model_number"]
    ax.ylabel = label_dict["mass"]
    lines!(ax, model_number * ones(length(mass)), mass, line_kwargs..., linewidth=1, color=burn)
    lines!(ax, model_number * ones(length(mass)), mass, line_kwargs..., linewidth=0.5, color=mixing)
end

"""
    function update_T_ρ_plot!(plot::StellarModels.JemsPlot, props::StellarModelProperties)

updates the observables of this Tρ `plot` with relevant data of the stellar model properties `props`.
"""
function update_Kipp_plot!(plot::StellarModels.JemsPlot, props::StellarModelProperties)
    make_Kipp_plot!(plot.ax, props.model_number, props.m[1:(props.nz)] / MSUN,
                    mixing_colors[get.(Ref(mixing_map), props.mixing_type[1:(props.nz)], missing)],
                    burning_map.(log10.(props.eps_nuc[1:(props.nz)])))
end
