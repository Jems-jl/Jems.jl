using Jems.NuclearNetworks

mutable struct AbundancePlot{TOBS, TOBS2, TAXIS}
    species_names::Vector{Symbol}
    xa_index::Dict{Symbol,Int}
    abundances::Vector{TOBS}
    label_locations::Vector{TOBS2}
    mass::TOBS
    axis::TAXIS
end

function AbundancePlot(grid_pos, network::TNET;
            log_yscale=false, ymin=0.0, ymax=1.1, colors=Makie.wong_colors(),
            label_fontsize=20) where{TNET<:AbstractNuclearNetwork}
    if !log_yscale
        axis = Axis(grid_pos, xlabel=L"\text{Mass}\;[M_\odot]", ylabel="Abundance")
    else
        axis = Axis(grid_pos, xlabel=L"\text{Mass}\;[M_\odot]", ylabel="Abundance", yscale=log10, yminorticks=IntervalsBetween(9))
    end

    mass = Observable(zeros(Float64,0))
    abundances = []
    label_locations = []
    for i in eachindex(network.species_names)
        color = colors[(i-1)%length(colors)+1]

        abundance_obs = Observable(zeros(Float64,0))
        push!(abundances, abundance_obs)

        lines!(axis, mass, abundance_obs, color=color)

        loc = Observable(Point(0.5,0.5))
        push!(label_locations, loc)

        text!(axis, loc, text=String(network.species_names[i]),
                color=color, align=(:left,:center),
                fontsize=label_fontsize)
    end

    ylims!(axis, ymin, ymax)

    # Here we do [abundances...] to get a typed vector. Same for label_locations
    return AbundancePlot(network.species_names, network.xa_index,
                            [abundances...], [label_locations...], mass, axis)
end

function update_plot!(p::AbundancePlot, m)
    p.mass[] = (@view m.props.m[1:m.props.nz]) / MSUN
    total_mass = m.props.m[m.props.nz]/MSUN
    for i in eachindex(p.species_names)
        index = p.xa_index[p.species_names[i]]
        p.abundances[i][] = [max(1e-99,get_value(m.props.xa[k, index])) for k in 1:m.props.nz]

        surf_value = get_value(m.props.xa[m.props.nz, index])
        p.label_locations[i][] = Point(total_mass, surf_value)
    end
    xlims!(p.axis, 0, total_mass*1.15)
end