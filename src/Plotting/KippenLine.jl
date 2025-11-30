mutable struct KippenLine{TOBS, TAXIS}
    xvals::TOBS
    yvals::TOBS
    xaxis::Symbol
    yaxis::Symbol
    xscale::Float64
    min_log_eps::Float64
    max_log_eps::Float64
    axis::TAXIS
end

function KippenLine(grid_pos; min_log_eps=0.0, max_log_eps=5.0, xaxis=:model_number, yaxis=:mass, time_units=:Myr)
    xscale=1.0
    if xaxis==:model_number
        xlabel="Model Number"
    elseif xaxis==:time
        if time_units==:yr
            xlabel=L"\text{Time}\;[\text{yr}]"
            xscale=SECYEAR
        elseif time_units==:kyr
            xlabel=L"\text{Time}\;[\text{kyr}]"
            xscale=SECYEAR*1e3
        elseif time_units==:Myr
            xlabel=L"\text{Time}\;[\text{Myr}]"
            xscale = SECYEAR*1e6
        elseif time_units==:Gyr
            xlabel=L"\text{Time}\;[\text{Gyr}]"
            xscale = SECYEAR*1e9
        else
            throw(ArgumentError("Invalid choice of time_units (:$time_units). Possible values are :yr, :kyr, :Myr, :Gyr"))
        end
    else
        throw(ArgumentError("Invalid choice of xaxis (:$xaxis). Possible values are :model_number, :time"))
    end
    if yaxis==:mass
        ylabel=L"\text{Mass}\;[M_\odot]"
    elseif yaxis==:radius
        ylabel=L"\text{Radius}\;[R_\odot]"
    else
        throw(ArgumentError("Invalid choice of yaxis (:$yaxis). Possible values are :mass, :radius"))
    end
    axis = Axis(grid_pos, xlabel=xlabel, ylabel=ylabel, ygridvisible=false, xgridvisible=false)
    xvals = Observable(zeros(Float64,0))
    yvals = Observable(zeros(Float64,0))
    lines!(axis, xvals, yvals, color=:black)
    return KippenLine(xvals, yvals, xaxis, yaxis, xscale, min_log_eps, max_log_eps, axis)
end

"""
    function draw_Kipp_lines!(p::KippenLine, model_number, mass, mixing, burn; line_kwargs=Dict())

Plots two lines per model in a Kippenhahn-like diagram (I dub it the KippenLine diagram).
One line contains info on the mixing state, the other on the burning regions.
"""
function draw_Kipp_lines!(p::KippenLine, xcoord, mass, mixing, burn; line_kwargs=Dict())
    lines!(p.axis, xcoord * ones(length(mass)), mass, line_kwargs..., linewidth=3, color=burn)
    # define a zvalue here as ones to ensure its plotted on top
    lines!(p.axis, xcoord * ones(length(mass)), mass, ones(length(mass)), line_kwargs..., linewidth=1, color=mixing)
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
function update_plot!(p::KippenLine, m)
    if p.xaxis == :model_number
        xcoord = m.props.model_number
    elseif p.xaxis == :time
        xcoord = m.props.time/p.xscale
    end
    p.xvals[] = push!(p.xvals[], xcoord)

    if p.yaxis == :mass
        p.yvals[] = push!(p.yvals[], m.props.mstar/MSUN)
        ycoords = (@view m.props.m[1:m.props.nz]) / MSUN
    elseif p.yaxis == :radius
        p.yvals[] = push!(p.yvals[], exp(get_value(m.props.lnr[m.props.nz]))/RSUN)
        ycoords = [exp(get_value(m.props.lnr[k]))/RSUN for k in 1:m.props.nz]
    end

    #draw the fun part
    draw_Kipp_lines!(p, xcoord, ycoords,
                     kipp_mixing_colors[(get.(Ref(mixing_map), (@view m.props.mixing_type[1:m.props.nz]), missing))],
                     burning_colors[burning_map.(log10.(abs.(@view m.props.ϵ_nuc[1:m.props.nz]));
                                                 min_log_eps=p.min_log_eps, max_log_eps=p.max_log_eps)])
end