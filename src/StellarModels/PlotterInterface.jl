using GLMakie

@kwdef mutable struct JemsPlot
    ax::Makie.Axis
    type::Symbol
    x_obs::Union{Dict{Symbol,Observable},Nothing}
    y_obs::Union{Dict{Symbol,Observable},Nothing}

    JemsPlot(ax, type) = new(ax, type, nothing, nothing)
end

mutable struct JemsFigure
    window::Makie.Figure
    plots::Vector{JemsPlot}
    scr::Union{MakieScreen,Nothing}

    JemsFigure(window, no_of_plots) = new(window, Vector{JemsPlot}(undef, no_of_plots), nothing)
        
end

"""
    mutable struct Plotter

Structure that contains references to windows, axes and observables to be plotted with Jems.Plotting
"""
@kwdef mutable struct Plotter
    figs::Union{Vector{JemsFigure},Nothing} = nothing
end
