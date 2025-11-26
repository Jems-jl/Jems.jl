export update_plotter!, AbstractPlotter

abstract type AbstractPlotter end

@kwdef struct Plotter{TFIG, TPLOTS} <: AbstractPlotter
    fig::TFIG
    plots::TPLOTS
    max_fps::Int
end

function update_plotter!(plotter::Plotter, m, refresh)
    for plot in plotter.plots
        update_plot!(plot, m)
    end
    if refresh
        display(plotter.fig)
        sleep(1/plotter.max_fps)
    end
end

struct NullPlotter <: AbstractPlotter end

function update_plotter!(plotter::NullPlotter, m, refresh)
end