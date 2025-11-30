mutable struct HRPlot{TOBS, TOBS2, TAXIS}
    logTeff::TOBS
    logL::TOBS
    current_position::TOBS2
    axis::TAXIS
end

function HRPlot(grid_pos)
    axis = Axis(grid_pos, xlabel=L"\log_{10}\left(T_\mathrm{eff}/[\mathrm{K}]\right)", ylabel=L"\log_{10}L/L_\odot", xreversed=true)
    logTeff = Observable(zeros(Float64,0))
    logL = Observable(zeros(Float64,0))
    lines!(axis, logTeff, logL)

    current_position = Observable(Point(0.0,0.0))
    scatter!(axis, current_position, marker=:circle, color=:red)

    return HRPlot(logTeff, logL, current_position, axis)
end

function update_plot!(p::HRPlot, m)
    new_logTeff = get_value(m.props.lnT[m.props.nz])*log10(ℯ)
    new_logL = log10(get_value(m.props.L[m.props.nz]))
    p.logTeff[] = push!(p.logTeff[], new_logTeff)
    p.logL[] = push!(p.logL[], new_logL)
    p.current_position[] = Point(new_logTeff, new_logL)
end