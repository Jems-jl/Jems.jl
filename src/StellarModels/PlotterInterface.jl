using GLMakie

"""
    mutable struct Plotter

Structure that contains references to windows, axes and observables to be plotted with Jems.Plotting
"""
@kwdef mutable struct Plotter
    fig::Union{Figure,Nothing} = nothing
    axs::Union{Vector{Makie.Block},Nothing} = nothing
    axno::Union{Dict{Symbol,Int},Nothing} = nothing
    obs::Union{Dict{Symbol, Union{Observable, Dict{Symbol, Observable}}},Nothing} = nothing
    scr::Union{MakieScreen,Nothing} = nothing
end

