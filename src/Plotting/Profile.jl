"""
    function create_profile_observables!(plot::StellarModels.JemsPlot, xvals::Dict, yvals::Dict;
                                         alt_yvals::Dict=nothing)

Creates relevant profile observables for this `plot` given `xvals`, `yvals` and 
optionally `alt_yvals` for the right hand axis.
"""
function create_profile_observables!(plot::StellarModels.JemsPlot, xvals::Dict, yvals::Dict; alt_yvals::Dict=nothing)
    for (name, vals) in pairs(xvals)
        plot.x_obs[name] = Observable{Vector{Float64}}(vals)
    end
    for (name, vals) in pairs(yvals)
        plot.y_obs[name] = Observable{Vector{Float64}}(vals)
    end
    if !isnothing(alt_yvals)
        for (name, vals) in pairs(alt_yvals)
            plot.alt_y_obs[name] = Observable{Vector{Float64}}(vals)
        end
    end
end

"""
    function make_profile_plot!(ax::Axis, xvals::Observable, yvals::Dict{Symbol,Observable};
                                xlabel::AbstractString="", ylabels::Dict{Symbol,<:AbstractString}=Dict(),
                                alt_ax::Axis=nothing, alt_yvals::Dict{Symbol,Observable}=nothing,
                                alt_ylabels::Dict{Symbol,<:AbstractString}=nothing,
                                line_kwargs=Dict())

Plots a line for each entry in the `yvals` dict, and puts the given `ylabels` in a legend.
The `alt_ax` is optional for the right hand side axis.
"""
function make_profile_plot!(ax::Axis, xvals::Observable, yvals::Dict{Symbol,Observable};
                            xlabel::AbstractString="", ylabels::Dict{Symbol,<:AbstractString}=Dict(),
                            alt_ax::Axis=nothing, alt_yvals::Dict{Symbol,Observable}=nothing,
                            alt_ylabels::Dict{Symbol,<:AbstractString}=nothing,
                            line_kwargs=Dict())
    ax.xlabel = xlabel
    (color, state) = iterate(colors)
    for (name, yline) in yvals
        lines!(ax, xvals, yline, line_kwargs..., color=color, label=ylabels[name])
        (color, state) = iterate(colors, state)
    end
    if !isnothing(alt_ax)
        for (name, yline) in alt_yvals
            lines!(alt_ax, xvals, yline, line_kwargs..., color=color, label=alt_ylabels[name])
            (color, state) = iterate(colors, state)
        end
        axislegend(ax, position=:lt)
        axislegend(alt_ax, position=:rt)
    else
        axislegend(ax)
    end
end

"""

    update_profile_plot!(plot::StellarModels.JemsPlot, sm::StellarModel)

updates the observables of this `plot` with relevant data of the stellar model `sm`.
"""
function update_profile_plot!(plot::StellarModels.JemsPlot, sm::StellarModel)
    for (key, obs) in pairs(plot.x_obs)
        obs.val = StellarModels.profile_output_functions[String(key)].((sm,), 1:(sm.props.nz))
    end
    for (key, obs) in pairs(plot.y_obs)
        obs.val = StellarModels.profile_output_functions[String(key)].((sm,), 1:(sm.props.nz))
    end
    for (key, obs) in pairs(plot.alt_y_obs)
        obs.val = StellarModels.profile_output_functions[String(key)].((sm,), 1:(sm.props.nz))
    end
end

