export update_plotter!, AbstractPlotter

abstract type AbstractPlotter end

@kwdef struct Plotter{TFIG, TPLOTS} <: AbstractPlotter
    fig::TFIG
    plots::TPLOTS
    max_fps::Int = 1000
    update_interval::Int = 1
    save_interval::Int = 0
    save_name::String = "plots_"
    save_folder::String = "png"
    model_number_pad::Int = 9
end

function update_plotter!(plotter::Plotter, m)
    for plot in plotter.plots
        update_plot!(plot, m)
    end
    if m.props.model_number % plotter.update_interval == 0
        display(plotter.fig)
        sleep(1/plotter.max_fps)
    end
    if plotter.save_interval > 0 &&  m.props.model_number % plotter.save_interval == 0
        if !isdir(plotter.save_folder)
            mkdir(plotter.save_folder)
        end
        save("$(plotter.save_folder)/$(plotter.save_name)$(string(m.props.model_number, pad=plotter.model_number_pad)).png", plotter.fig)
    end
end

struct NullPlotter <: AbstractPlotter end

function update_plotter!(plotter::NullPlotter, m, refresh)
end