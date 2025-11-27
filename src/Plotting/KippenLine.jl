mutable struct KippenLine{TOBS, TAXIS}
    xvals::TOBS
    yvals::TOBS
    xaxis::Symbol
    yaxis::Symbol
    xaxis_scale::Float64
    min_log_eps::Float64
    max_log_eps::Float64
    axis::TAXIS
end

function KippenLine(grid_pos; min_log_eps=0.0, max_log_eps=5.0, xaxis=:model_number, yaxis=:mass, xaxis_scale=1.0)
    if xaxis==:model_number
        xlabel="Model Number"
    elseif xaxis==:time
        xlabel="Time"
    end
    if yaxis==:mass
        ylabel="Mass"
    elseif yaxis==:radius
        ylabel="Radius"
    end
    axis = Axis(grid_pos, xlabel=xlabel, ylabel=ylabel, ygridvisible=false, xgridvisible=false)
    xvals_start::Vector{Float64} = []
    yvals_start::Vector{Float64} = []
    xvals = Observable(xvals_start)
    yvals = Observable(yvals_start)
    lines!(axis, xvals, yvals, color=:black)
    return KippenLine(xvals, yvals, xaxis, yaxis, xaxis_scale, min_log_eps, max_log_eps, axis)
end

"""
    function draw_Kipp_lines!(p::KippenLine, model_number, mass, mixing, burn; line_kwargs=Dict())

Plots two lines per model in a Kippenhahn-like diagram (I dub it the KippenLine diagram).
One line contains info on the mixing state, the other on the burning regions.
"""
function draw_Kipp_lines!(p::KippenLine, model_number, mass, mixing, burn; line_kwargs=Dict())
    lines!(p.axis, model_number * ones(length(mass)), mass, line_kwargs..., linewidth=1, color=burn)
    lines!(p.axis, (model_number + 0.5) * ones(length(mass)), mass, line_kwargs..., linewidth=1, color=mixing)
end

function burning_map(log_eps_nuc; min_log_eps=0.0, max_log_eps=15.0)  # map log eps nuc to interval [0.0, 1.0]
    if log_eps_nuc < min_log_eps
        return 0.0
    elseif log_eps_nuc > max_log_eps
        return 1.0
    else
        return log_eps_nuc / (max_log_eps - min_log_eps)
    end
end

kipp_mixing_colors = [RGBAf(0.5, 0.5, 0.5), RGBAf(0, 0, 1)]
kipp_mixing_colors[1] = RGBAf(1, 1, 1)  # make no_mixing white in KippenLine diagram to avoid clutter
burning_colors = cgrad(:linear_wyor_100_45_c55_n256)
const mixing_map = Dict(:no_mixing => 1,
                        :convection => 2)
function update_plot!(p::KippenLine, m)
    if p.xaxis == :model_number
        p.xvals[] = push!(p.xvals[], m.props.model_number)
    # do time
    end
    if p.yaxis == :mass
        p.yvals[] = push!(p.yvals[], m.props.mstar/MSUN)
    # do radius
    end

    #draw the fun part
    draw_Kipp_lines!(p, m.props.model_number, (@view m.props.m[1:m.props.nz]) / MSUN,
                     kipp_mixing_colors[(get.(Ref(mixing_map), (@view m.props.mixing_type[1:m.props.nz]), missing))],
                     burning_colors[burning_map.(log10.(abs.(@view m.props.ϵ_nuc[1:m.props.nz]));
                                                 min_log_eps=p.min_log_eps, max_log_eps=p.max_log_eps)])
end