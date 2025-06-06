using GLMakie

@kwdef mutable struct JemsPlot
    ax::Makie.Axis
    type::Symbol
    x_obs::Union{Dict{Symbol,Observable},Nothing}
    y_obs::Union{Dict{Symbol,Observable},Nothing}
    other_obs::Union{Dict{Symbol,Observable},Nothing}
    # for right-hand-side axes
    alt_ax::Union{Makie.Axis,Nothing}
    alt_y_obs::Union{Dict{Symbol,Observable},Nothing}

    JemsPlot(ax, type) = new(ax, type, nothing, nothing, nothing, nothing, nothing)
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
