"""
    init_plots(sm::StellarModel)

Sets up all observables to be traced this run, and creates the figure and axis where they will be plotted
"""
function init_plots!(sm::StellarModel)
    basic_theme = Theme(fonts=(regular=texfont(:text), bold=texfont(:bold),
                               italic=texfont(:italic), bold_italic=texfont(:bolditalic)),
                        fontsize=50, size=(1000, 750), linewidth=7,
                        Axis=(xlabelsize=40, ylabelsize=40, titlesize=40, xgridvisible=false, ygridvisible=false,
                              spinewidth=2.5, xminorticksvisible=true, yminorticksvisible=true, xtickalign=1,
                              ytickalign=1,
                              xminortickalign=1, yminortickalign=1, xticksize=14, xtickwidth=2.5, yticksize=14,
                              ytickwidth=2.5, xminorticksize=7, xminortickwidth=2.5, yminorticksize=7,
                              yminortickwidth=2.5,
                              xticklabelsize=35, yticklabelsize=35, xticksmirrored=true, yticksmirrored=true),
                        Legend=(patchsize=(70, 10), framevisible=false, patchlabelgap=20, rowgap=10))

    GLMakie.set_theme!(basic_theme)
    GLMakie.activate!()
    #GLMakie.set_window_config!(; float=true)  # place windows on top # this does not work in newer GLMakie versions it seems

    # create figure/axes objects
    init_figure!(sm)

    # populate observables and plot elements on the respective axes
    for plot in sm.plt.plots
        plot.x_obs = Dict{Symbol,Observable}()
        plot.y_obs = Dict{Symbol,Observable}()
        plot.alt_y_obs = Dict{Symbol,Observable}()
        if plot.type == :HR
            create_HR_observables!(plot, sm.props)
            make_HR_plot!(plot.ax, plot.x_obs[:Teff], plot.y_obs[:L], plot.x_obs[:Teff_now],
                          plot.y_obs[:L_now]; scatter_kwargs=Dict(:color => "red", :markersize => 20))
        elseif plot.type == :profile
            # make observables
            xname = sm.opt.plotting.profile_xaxis
            xvals = Dict(Symbol(xname) => StellarModels.profile_output_functions[xname].((sm,), 1:(sm.nz)))
            ynames = sm.opt.plotting.profile_yaxes
            yvals = Dict([Symbol(name) =>
                        StellarModels.profile_output_functions[name].((sm,), 1:(sm.nz)) for name in ynames])
            altynames = sm.opt.plotting.profile_alt_yaxes
            if !isnothing(plot.alt_ax)
                altyvals = Dict([Symbol(name) =>
                            StellarModels.profile_output_functions[name].((sm,), 1:(sm.nz)) for name in altynames])
            else
                altyvals = nothing
            end
            create_profile_observables!(plot, xvals, yvals; alt_yvals=altyvals)
            # set labels
            ylabels = Dict([Symbol(name) => label_dict[name] for name in ynames])
            if !isnothing(plot.alt_ax)
                alt_ylabels = Dict([Symbol(name) => label_dict[name] for name in altynames])
            else
                alt_ylabels = nothing
            end
            make_profile_plot!(plot.ax, collect(values(plot.x_obs))[1], plot.y_obs;
                                xlabel=label_dict[sm.opt.plotting.profile_xaxis], ylabels=ylabels,
                                alt_ax=plot.alt_ax, alt_yvals=plot.alt_y_obs, alt_ylabels=alt_ylabels)
        elseif plot.type == :history
            xname = sm.opt.plotting.history_xaxis
            xvals = Dict(Symbol(xname) => StellarModels.history_output_functions[xname](sm))
            
            ynames = sm.opt.plotting.history_yaxes
            yvals = Dict([Symbol(name) => StellarModels.history_output_functions[name](sm) for name in ynames])
            
            altynames = sm.opt.plotting.history_alt_yaxes
            if !isnothing(plot.alt_ax)
                altyvals = Dict([Symbol(name) => StellarModels.history_output_functions[name](sm) for name in altynames])
            else
                altyvals = nothing
            end
            create_history_observables!(plot, xvals, yvals; alt_yvals=altyvals)
            # set labels
            ylabels = Dict([Symbol(name) => label_dict[name] for name in ynames])
            if !isnothing(plot.alt_ax)
                alt_ylabels = Dict([Symbol(name) => label_dict[name] for name in altynames])
            else
                alt_ylabels = nothing
            end
            make_history_plot!(plot.ax, collect(values(plot.x_obs))[1], plot.y_obs;
                                xlabel=label_dict[sm.opt.plotting.history_xaxis], ylabels=ylabels,
                                alt_ax=plot.alt_ax, alt_yvals=plot.alt_y_obs, alt_ylabels=alt_ylabels)
        end
    end

    # display the figure
    sm.plt.scr = display(sm.plt.fig)
end

"""
    init_figure!(sm::StellarModel)

Initializes the plotting figure and adds axes
"""
function init_figure!(sm::StellarModel)
    sm.plt.fig = Figure()
    no_of_plots = length(sm.opt.plotting.window_specs)
    sm.plt.plots = Vector{StellarModels.JemsPlot}(undef, no_of_plots)
    for j = 1:no_of_plots
        this_axis = Axis(sm.plt.fig[sm.opt.plotting.window_layouts[j]...])
        this_type = Symbol(sm.opt.plotting.window_specs[j])
        sm.plt.plots[j] = StellarModels.JemsPlot(this_axis, this_type)
        
        if (sm.plt.plots[j].type == :profile && length(sm.opt.plotting.profile_alt_yaxes) > 0) || 
            (sm.plt.plots[j].type == :history && length(sm.opt.plotting.history_alt_yaxes) > 0)
            # modifications if we have an alt axis
            sm.plt.plots[j].alt_ax = Axis(sm.plt.fig[sm.opt.plotting.window_layouts[j]...], yaxisposition = :right)
            hidespines!(sm.plt.plots[j].alt_ax, :l, :t, :b)
            hidexdecorations!(sm.plt.plots[j].alt_ax)
            hidespines!(sm.plt.plots[j].ax, :r)
            sm.plt.plots[j].ax.yticksmirrored = false
            sm.plt.plots[j].alt_ax.yticksmirrored = false
            linkxaxes!(sm.plt.plots[j].ax, sm.plt.plots[j].alt_ax)
        end
    end
end