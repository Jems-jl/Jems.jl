using Jems.StellarModels # label names and functions come from IO.jl in StellarModels

mutable struct ProfilePlot{TOBS, TAXIS, TAXIS2}
    xvals::TOBS
    yvals::TOBS
    other_yvals::TOBS
    function_x::Function
    function_y::Function
    function_othery::Union{Function,Nothing}
    axis::TAXIS
    other_axis::TAXIS2
end

function ProfilePlot(grid_pos, sm; x_name="", y_name="", othery_name="", link_yaxes=false, ycolor=nothing, othery_color=nothing)
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
    xlabel = sm.profile_output_labels[x_name]
    ylabel = sm.profile_output_labels[y_name]
    axis = Axis(grid_pos, xlabel=xlabel, ylabel=ylabel, xgridvisible=false, ygridvisible=false, ylabelcolor=ylabelcolor,
                    yticksmirrored=yticksmirrored)
    if othery_name != ""
        othery_label = sm.profile_output_labels[othery_name]
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
    function_x = sm.profile_output_functions[x_name]
    function_y = sm.profile_output_functions[y_name]
    if othery_name != ""
        function_othery = sm.profile_output_functions[othery_name]
    else
        function_othery = nothing
    end

    return ProfilePlot(xvals, yvals, other_yvals,
                        function_x, function_y, function_othery, axis, other_axis)
end

function update_plot!(p::ProfilePlot, m)
    p.xvals[] = [p.function_x(m,k) for k in 1:m.props.nz]
    p.yvals[] = [p.function_y(m,k) for k in 1:m.props.nz]
    if !isnothing(p.other_axis)
        p.other_yvals[] = [Float64(p.function_othery(m,k)) for k in 1:m.props.nz]
    end
end