using GLMakie

@kwdef mutable struct Plotter
    fig::Union{Figure,Nothing}
    axs::Union{Vector{Makie.Block},Nothing}
    obs::Union{Dict{Symbol,Observable},Nothing}
    scr::Union{MakieScreen,Nothing}
end

function Plotter(x::Nothing)
    Plotter(fig=x, axs=x, obs=x, scr=x)
end
