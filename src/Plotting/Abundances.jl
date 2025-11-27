using Jems.NuclearNetworks

mutable struct AbundancePlot{TOBS, TAXIS}
    species_names::Vector{Symbol}
    xa_index::Dict{Symbol,Int}
    abundances::Vector{TOBS}
    mass::TOBS
    axis::TAXIS
end

function AbundancePlot(grid_pos, network::TNET; log_yscale=false, ymin=0.0, ymax=1.0) where{TNET<:AbstractNuclearNetwork}
    if !log_yscale
        axis = Axis(grid_pos, xlabel="Mass", ylabel="abundance")
    else
        axis = Axis(grid_pos, xlabel="Mass", ylabel="abundance", yscale=log10, yminorticks=IntervalsBetween(9))
    end

    mass = Observable(zeros(Float64,0))
    abundances = []
    for i in eachindex(network.species_names)
        abundance_obs = Observable(zeros(Float64,0))
        push!(abundances, abundance_obs)

        lines!(axis, mass, abundance_obs)
    end

    ylims!(axis, ymin, ymax)

    # Here we do [abundances...] to get a typed vector
    return AbundancePlot(network.species_names, network.xa_index, [abundances...], mass, axis)
end

function update_plot!(p::AbundancePlot, m)
    p.mass[] = (@view m.props.m[1:m.props.nz]) / MSUN
    for i in eachindex(p.species_names)
        index = p.xa_index[p.species_names[i]]
        p.abundances[i][] = [max(1e-99,get_value(m.props.xa[k, index])) for k in 1:m.props.nz]
    end
end