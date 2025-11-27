mutable struct HRPlot{TOBS, TAXIS}
    logTeff::TOBS
    logL::TOBS
    axis::TAXIS
end

function HRPlot(grid_pos)
    axis = Axis(grid_pos, xlabel="logTeff", ylabel="logL")
    logTeff_start::Vector{Float64} = []
    logL_start::Vector{Float64} = []
    logTeff = Observable(logTeff_start)
    logL = Observable(logL_start)
    lines!(axis, logTeff, logL)
    return HRPlot(logTeff, logL, axis)
end

function update_plot!(p::HRPlot, m)
    new_logTeff = get_value(m.props.lnT[m.props.nz])*log10(ℯ)
    new_logL = log10(get_value(m.props.L[m.props.nz]))
    p.logTeff[] = push!(p.logTeff[], new_logTeff)
    p.logL[] = push!(p.logL[], new_logL)
end