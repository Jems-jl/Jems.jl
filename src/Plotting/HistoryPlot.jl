using Jems.StellarModels # label names and functions come from IO.jl in StellarModels

mutable struct HistoryPlot{TOBS, TAXIS, TAXIS2}
    xvals::TOBS
    yvals::TOBS
    other_yvals::TOBS
    function_x::Function
    function_y::Function
    function_othery::Union{Function,Nothing}
    axis::TAXIS
    other_axis::TAXIS2
end

function HistoryPlot(grid_pos, sm::StellarModel; x_name="", y_name="", othery_name="", link_yaxes=false, ycolor=nothing, othery_color=nothing)
    if isnothing(ycolor)
        ycolor = Makie.wong_colors()[1]
    end
    if isnothing(othery_color)
        othery_color = Makie.wong_colors()[2]
    end
    # TODO check if key is available and provide instructive error message if not
    if othery_name != ""
        yticksmirrored = false
        ylabelcolor = ycolor
    else
        yticksmirrored = true
        ylabelcolor = :black
    end
    xlabel = sm.history_output_labels[x_name]
    ylabel = sm.history_output_labels[y_name]
    axis = Axis(grid_pos, xlabel=xlabel, ylabel=ylabel, xgridvisible=false, ygridvisible=false, ylabelcolor=ylabelcolor,
                    yticksmirrored=yticksmirrored)
    if othery_name != ""
        othery_label = sm.history_output_labels[othery_name]
        other_axis = Axis(grid_pos, yaxisposition = :right, ylabel=othery_label, ylabelcolor = othery_color,
                            xgridvisible=false, ygridvisible=false, yticksmirrored=yticksmirrored)
        hidespines!(other_axis)
        hidexdecorations!(other_axis)
        if link_yaxes
            linkyaxes!(axis, other_axis)
        end
    else
        other_axis = nothing
    end

    xvals = Observable(zeros(Float64,0))
    yvals = Observable(zeros(Float64,0))
    other_yvals = Observable(zeros(Float64,0))
    lines!(axis, xvals, yvals, color=ycolor)
    if othery_name != ""
        lines!(other_axis, xvals, other_yvals, color=othery_color)
    end
    function_x = sm.history_output_functions[x_name]
    function_y = sm.history_output_functions[y_name]
    if othery_name != ""
        function_othery = sm.history_output_functions[othery_name]
    else
        function_othery = nothing
    end

    return HistoryPlot(xvals, yvals, other_yvals,
                        function_x, function_y, function_othery, axis, other_axis)
end

function update_plot!(p::HistoryPlot, m)
    p.xvals[] = push!(p.xvals[], p.function_x(m))
    p.yvals[] = push!(p.yvals[], p.function_y(m))
    if !isnothing(p.other_axis)
        p.other_yvals[] = push!(p.other_yvals[], p.function_othery(m))
    end
end