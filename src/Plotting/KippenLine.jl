"""
    function create_Kipp_observables!(plot::StellarModels.JemsPlot, props::StellarModelProperties)

Creates relevant observables for this `plot` given the StellarModelProperties `props`
"""
function create_Kipp_observables!(plot::StellarModels.JemsPlot, props::StellarModelProperties)
    plot.y_obs[:mass] = Observable(Float64[])
    push!(plot.y_obs[:mass][], props.mstar / MSUN)
    plot.x_obs[:model_number] = Observable(Int64[])
    push!(plot.x_obs[:model_number][], props.model_number)
end

"""
    function draw_Kipp_lines!(ax::Axis)

Plots two lines per model in a Kippenhahn-like diagram (I dub it the KippenLine diagram).
One line contains info on the mixing state, the other on the burning regions.
"""
function draw_Kipp_lines!(ax::Axis, model_number, mass, mixing, burn; line_kwargs=Dict())
    lines!(ax, model_number * ones(length(mass)), mass, line_kwargs..., linewidth=1, color=burn)
    lines!(ax, (model_number + 0.5) * ones(length(mass)), mass, line_kwargs..., linewidth=1, color=mixing)
end

function make_Kipp_plot!(ax::Axis, model_obs::Observable, mass_obs::Observable, line_kwargs=Dict())
    ax.xlabel = label_dict["model_number"]
    ax.ylabel = label_dict["mass"]
    lines!(ax, model_obs, mass_obs, line_kwargs..., color=:gray)
    translate!(Accum, ax.scene.plots[1], Vec3f(0, 0, 100))  # raise this line so it renders above the KippLines
end

"""
    function update_T_ρ_plot!(plot::StellarModels.JemsPlot, props::StellarModelProperties)

updates the observables of this Tρ `plot` with relevant data of the stellar model properties `props`.
"""
function update_Kipp_plot!(plot::StellarModels.JemsPlot, props::StellarModelProperties,
                           plotopt::StellarModels.PlottingOptions)
    draw_Kipp_lines!(plot.ax, props.model_number, (@view props.m[1:(props.nz)]) / MSUN,
                     kipp_mixing_colors[(get.(Ref(mixing_map), (@view props.mixing_type[1:(props.nz)]), missing))],
                     burning_colors[burning_map.(log10.(abs.(@view props.eps_nuc[1:(props.nz)]));
                                                 min_log_eps=plotopt.min_log_eps, max_log_eps=plotopt.max_log_eps)])
    push!(plot.y_obs[:mass][], props.mstar / MSUN)
    push!(plot.x_obs[:model_number][], props.model_number)
end
