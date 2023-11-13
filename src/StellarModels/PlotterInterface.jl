using GLMakie

@kwdef mutable struct JemsPlot
    ax::Makie.Axis
    type::Symbol
    x_obs::Union{Dict{Symbol,Observable},Nothing}
    y_obs::Union{Dict{Symbol,Observable},Nothing}

    JemsPlot(ax, type) = new(ax, type, nothing, nothing)
end

"""
    mutable struct Plotter

Structure that contains references to windows, axes and observables to be plotted with Jems.Plotting
"""
@kwdef mutable struct Plotter
    fig::Union{Makie.Figure, Nothing} = nothing
    plots::Union{Vector{JemsPlot}, Nothing} = nothing
    scr::Union{MakieScreen,Nothing} = nothing
end
