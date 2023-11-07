using GLMakie

"""
    mutable struct Plotter

Structure that contains references to windows, axes and observables to be plotted with Jems.Plotting
"""
@kwdef mutable struct Plotter
    fig::Union{Figure,Nothing}
    axs::Union{Vector{Makie.Block},Nothing}
    obs::Union{Dict{Symbol,Observable},Nothing}
    scr::Union{MakieScreen,Nothing}
end

"""
    Plotter constructor
"""
function Plotter(x::Nothing)
    Plotter(fig=x, axs=x, obs=x, scr=x)
end
