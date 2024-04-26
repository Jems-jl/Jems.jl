"""
    init_plots(sm::StellarModel)

Sets up all observables to be traced this run, and creates the figure and axis where they will be plotted,
executes after the first model is found
"""
function init_plots!(m::AbstractModel)
    basic_theme = Theme(fonts=(regular=texfont(:text), bold=texfont(:bold),
                               italic=texfont(:italic), bold_italic=texfont(:bolditalic)),
                        fontsize=30, size=(1000, 750), linewidth=7,
                        Axis=(xlabelsize=40, ylabelsize=40, titlesize=40, xgridvisible=false, ygridvisible=false,
                              spinewidth=2.5, xminorticksvisible=true, yminorticksvisible=true, xtickalign=1,
                              ytickalign=1,
                              xminortickalign=1, yminortickalign=1, xticksize=14, xtickwidth=2.5, yticksize=14,
                              ytickwidth=2.5, xminorticksize=7, xminortickwidth=2.5, yminorticksize=7,
                              yminortickwidth=2.5,
                              xticklabelsize=35, yticklabelsize=35, xticksmirrored=true, yticksmirrored=true),
                        Legend=(patchsize=(70, 10), framevisible=false, patchlabelgap=20, rowgap=10, fontsize=12))

    GLMakie.set_theme!(basic_theme)
    GLMakie.activate!(fxaa=false, ssao=false)
    #GLMakie.set_window_config!(; float=true)  # place windows on top # this does not work in newer GLMakie versions it seems

    # create figure/axes objects
    init_figure!(m)

    # populate observables and plot elements on the respective axes
    for plot in m.plt.plots
        plot.x_obs = Dict{Symbol,Observable}()
        plot.y_obs = Dict{Symbol,Observable}()
        if plot.type == :HR
            create_HR_observables!(plot, m.props)
            make_HR_plot!(plot.ax, plot.x_obs[:Teff], plot.y_obs[:L], plot.x_obs[:Teff_now],
                          plot.y_obs[:L_now]; scatter_kwargs=Dict(:color => "red", :markersize => 20))
        elseif plot.type == :profile
            # make observables
            xname = m.opt.plotting.profile_xaxis
            xvals = Dict(Symbol(xname) => StellarModels.profile_output_functions[xname].((m,), 1:(m.props.nz)))
            ynames = m.opt.plotting.profile_yaxes
            yvals = Dict([Symbol(name) => StellarModels.profile_output_functions[name].((m,), 1:(m.props.nz))
                          for name in ynames])
            altynames = m.opt.plotting.profile_alt_yaxes
            if !isnothing(plot.alt_ax)
                plot.alt_y_obs = Dict{Symbol,Observable}()
                altyvals = Dict([Symbol(name) => StellarModels.profile_output_functions[name].((m,), 1:(m.props.nz))
                                 for name in altynames])
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
                               xlabel=label_dict[m.opt.plotting.profile_xaxis], ylabels=ylabels,
                               alt_ax=plot.alt_ax, alt_yvals=plot.alt_y_obs, alt_ylabels=alt_ylabels)
        elseif plot.type == :TRhoProfile
            plot.other_obs = Dict{Symbol,Observable}()
            create_T_ρ_observables!(plot, m.props)
            make_T_ρ_plot!(plot.ax, plot.x_obs[:log_ρ], plot.y_obs[:log_T], plot.other_obs[:colors])
        elseif plot.type == :history
            xname = m.opt.plotting.history_xaxis
            xvals = Dict(Symbol(xname) => StellarModels.history_output_functions[xname](m))

            ynames = m.opt.plotting.history_yaxes
            yvals = Dict([Symbol(name) => StellarModels.history_output_functions[name](m) for name in ynames])

            altynames = m.opt.plotting.history_alt_yaxes
            if !isnothing(plot.alt_ax)
                plot.alt_y_obs = Dict{Symbol,Observable}()
                altyvals = Dict([Symbol(name) => StellarModels.history_output_functions[name](m) for name in altynames])
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
                               xlabel=label_dict[m.opt.plotting.history_xaxis], ylabels=ylabels,
                               alt_ax=plot.alt_ax, alt_yvals=plot.alt_y_obs, alt_ylabels=alt_ylabels)
        elseif plot.type == :Kippenhahn
            create_Kipp_observables!(plot, m.props)
            make_Kipp_plot!(plot.ax, plot.x_obs[:model_number], plot.y_obs[:mass])
            draw_Kipp_lines!(plot.ax, m.props.model_number, (@view m.props.m[1:(m.props.nz)]) / MSUN,
                             kipp_mixing_colors[(get.(Ref(mixing_map), (@view m.props.mixing_type[1:(m.props.nz)]),
                                                      missing))],
                             burning_colors[burning_map.(log10.(abs.(@view m.props.ϵ_nuc[1:(m.props.nz)]));
                                                         min_log_eps=m.opt.plotting.min_log_eps,
                                                         max_log_eps=m.opt.plotting.max_log_eps)])
        end
    end

    # display the figure
    m.plt.scr = display(m.plt.fig)
end

"""
    init_figure!(m<:AbstractModel)

Initializes the plotting figure and adds axes
"""
function init_figure!(m::AbstractModel)
    m.plt.fig = Figure()
    no_of_plots = length(m.opt.plotting.window_specs)
    m.plt.plots = Vector{StellarModels.JemsPlot}(undef, no_of_plots)
    for j = 1:no_of_plots
        if j <= length(m.opt.plotting.yaxes_log) && m.opt.plotting.yaxes_log[j]
            scale = log10
        else
            scale = identity
        end
        this_axis = Axis(m.plt.fig[m.opt.plotting.window_layout[j]...], yscale=scale)
        this_type = Symbol(m.opt.plotting.window_specs[j])
        m.plt.plots[j] = StellarModels.JemsPlot(this_axis, this_type)

        if (m.plt.plots[j].type == :profile && length(m.opt.plotting.profile_alt_yaxes) > 0) ||
           (m.plt.plots[j].type == :history && length(m.opt.plotting.history_alt_yaxes) > 0)
            # modifications if we have an alt axis
            if j <= length(m.opt.plotting.alt_yaxes_log) && m.opt.plotting.alt_yaxes_log[j]
                scale = log10
            else
                scale = identity
            end
            m.plt.plots[j].alt_ax = Axis(m.plt.fig[m.opt.plotting.window_layout[j]...], yaxisposition=:right,
                                         yscale=scale)
            hidespines!(m.plt.plots[j].alt_ax, :l, :t, :b)
            hidexdecorations!(m.plt.plots[j].alt_ax)
            hidespines!(m.plt.plots[j].ax, :r)
            m.plt.plots[j].ax.yticksmirrored = false
            m.plt.plots[j].alt_ax.yticksmirrored = false
            linkxaxes!(m.plt.plots[j].ax, m.plt.plots[j].alt_ax)
        end
    end
end