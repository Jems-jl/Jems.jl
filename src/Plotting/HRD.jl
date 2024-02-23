using Jems.DualSupport

"""
    create_HR_observables!(sm::StellarModel, plot::StellarModels.JemsPlot)

Creates teff and L observables and adds them to the observable list of the given plot
"""
function create_HR_observables!(plot::StellarModels.JemsPlot, props::StellarModelProperties)
    teff = exp(get_cell_value(props.lnT[props.nz]))
    plot.x_obs[:Teff_now] = Observable{Float64}(teff)
    plot.x_obs[:Teff] = Observable(Float64[])
    push!(plot.x_obs[:Teff][], teff)
    plot.y_obs[:L_now] = Observable{Float64}(get_cell_value(props.L[props.nz]))
    plot.y_obs[:L] = Observable(Float64[])
    push!(plot.y_obs[:L][], get_cell_value(props.L[props.nz]))
end

"""
    init_HR_plot!(ax::Axis, Teff::Observable, L::Observable, Teff_now::Observable, L_now::Observable;
                      line_kwargs=Dict(), scatter_kwargs=Dict())

Sets up the plot elements for an HRD
"""
function make_HR_plot!(ax::Axis, Teff::Observable, L::Observable, Teff_now::Observable, L_now::Observable;
                       line_kwargs=Dict(), scatter_kwargs=Dict())
    ax.xlabel = L"T_\mathrm{eff}/K"
    ax.ylabel = L"L/L_\odot"
    ax.xreversed = true
    lines!(ax, Teff, L; line_kwargs...)
    scatter!(ax, Teff_now, L_now; scatter_kwargs...)
end

"""
    function update_HR_plot!(plot::StellarModels.JemsPlot, sm::StellarModel)

Updates the given `plot` with the relevant HR data from the properties of the stellar model `props`.
"""
function update_HR_plot!(plot::StellarModels.JemsPlot, props::StellarModelProperties)
    push!(plot.x_obs[:Teff].val, exp(get_cell_value(props.lnT[props.nz])))
    plot.x_obs[:Teff_now].val = exp(get_cell_value(props.lnT[props.nz]))
    push!(plot.y_obs[:L].val, get_cell_value(props.L[props.nz]))
    plot.y_obs[:L_now].val = get_cell_value(props.L[props.nz])
end
