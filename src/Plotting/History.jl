"""
    create_history_observables!(sm::StellarModel, plot::StellarModels.JemsPlot)

creates the x and y observables and adds them to the observable list of the given plot
"""
function create_history_observables!(plot::StellarModels.JemsPlot, xval::Dict, yvals::Dict; alt_yvals::Dict=nothing)
    xkey = collect(keys(xval))[1]
    plot.x_obs[xkey] = Observable{Vector{Float64}}([0.0])
    plot.x_obs[xkey][][1] = xval[xkey]
    for (name, val) in pairs(yvals)
        plot.y_obs[name] = Observable{Vector{Float64}}([0.0])
        plot.y_obs[name][][1] = val
    end
    if !isnothing(alt_yvals)
        for (name, val) in pairs(alt_yvals)
            plot.alt_y_obs[name] = Observable{Vector{Float64}}([0.0])
            plot.alt_y_obs[name][][1] = val
        end
    end 
end

"""
    init_HR_plot!(ax::Axis, xval::Observable, yval::Observable; line_kwargs=Dict())

Sets up the plot elements for an history plot
"""
function make_history_plot!(ax::Axis, xval::Observable, yvals::Dict{Symbol,Observable};
                            xlabel::AbstractString="", ylabels::Dict{Symbol,<:AbstractString}=Dict(),
                            alt_ax::Axis=nothing, alt_yvals::Dict{Symbol,Observable}=nothing,
                            alt_ylabels::Dict{Symbol,<:AbstractString}=nothing,
                            line_kwargs=Dict())
    ax.xlabel = xlabel
    (color, state) = iterate(colors)
    for (name, yval) in pairs(yvals)
        lines!(ax, xval, yval; line_kwargs..., color=color, label=ylabels[name])
        (color, state) = iterate(colors, state)
    end
    if !isnothing(alt_ax)
        for (name, alt_yval) in pairs(alt_yvals)
            lines!(alt_ax, xval, alt_yval; line_kwargs..., color=color, label=alt_ylabels[name])
            (color, state) = iterate(colors, state)
        end
        axislegend(ax, position=:lt)
        axislegend(alt_ax, position=:rt)
    else
        axislegend(ax)
    end
end

function update_history_plot!(plot::StellarModels.JemsPlot, sm::StellarModel)
    for (key, obs) in pairs(plot.x_obs)
        push!(obs.val, StellarModels.history_output_options[String(key)][2](sm))
    end
    for (key, obs) in pairs(plot.y_obs)
        push!(obs.val, StellarModels.history_output_options[String(key)][2](sm))
    end
    for (key, obs) in pairs(plot.alt_y_obs)
        push!(obs.val, StellarModels.history_output_options[String(key)][2](sm))
    end
end
