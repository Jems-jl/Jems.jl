export update_plotter!, AbstractPlotter

abstract type AbstractPlotter end

@kwdef struct Plotter{TFIG, TPLOTS} <: AbstractPlotter
    fig::TFIG
    plots::TPLOTS
    max_fps::Int
    update_interval::Int
end

function update_plotter!(plotter::Plotter, m)
    for plot in plotter.plots
        update_plot!(plot, m)
    end
    if m.props.model_number % plotter.update_interval == 0
        display(plotter.fig)
        sleep(1/plotter.max_fps)
    end
end

struct NullPlotter <: AbstractPlotter end

function update_plotter!(plotter::NullPlotter, m, refresh)
end