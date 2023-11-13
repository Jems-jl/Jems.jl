module Plotting

using GLMakie, LaTeXStrings, MathTeXEngine, Jems.StellarModels
const colors = Iterators.cycle([:red, :blue, :green])
const label_dict = Dict("mass" => L"m / M_\odot",
                        "X" => L"X",
                        "Y" => L"Y",
                        "log10_T" => L"\log_{10}(T / K)",
                        "zone" => L"\mathrm{zone}")

"""
    init_plots(sm::StellarModel)

Sets up all observables to be traced this run, and creates the figure and axis where they will be plotted
"""
function init_plots!(sm::StellarModel)
    basic_theme = Theme(fonts=(regular=texfont(:text), bold=texfont(:bold),
                               italic=texfont(:italic), bold_italic=texfont(:bolditalic)),
                        fontsize=30, resolution=(1000, 750), linewidth=7,
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
    GLMakie.set_window_config!(; float=true)  # place windows on top

    # create figure/axes objects
    init_figures!(sm)

    # populate observables and plot elements on the respective axes
    for plot in sm.plt.plots
        if plot.type == :HR
            create_HR_observables!(sm, plot)
            make_HR_plot!(plot.ax, plot.x_obs[:Teff], plot.y_obs[:L], plot.x_obs[:Teff_now],
                        plot.y_obs[:L_now]; scatter_kwargs=Dict(:color => "red", :markersize => 20))
        elseif plot.type == :profile
            create_profile_observables!(sm, plot)
            ylabels = Dict{Symbol,String}()
            for name in sm.opt.plotting.profile_yaxes
                ylabels[Symbol(name)] = label_dict[name]
            end
            if !isnothing(plot.alt_ax)
                alt_ylabels = Dict{Symbol,String}()
                for name in sm.opt.plotting.profile_alt_yaxes
                    alt_ylabels[Symbol(name)] = label_dict[name]
                end
            else
                alt_ylabels = nothing
            end
            make_profile_plot!(plot.ax, plot.x_obs[:profile_x], plot.y_obs, label_dict[sm.opt.plotting.profile_xaxis], 
                            ylabels; alt_ax=plot.alt_ax, alt_yvals=plot.alt_y_obs, alt_ylabels=alt_ylabels )
        end
    end
    
    # display the figures
    sm.plt.scr = display(sm.plt.fig)
end

"""
    init_figure!(sm::StellarModel)

Initializes the plotting figures and adds axes
"""
function init_figures!(sm::StellarModel)
    sm.plt.fig = Figure()
    no_of_plots = length(sm.opt.plotting.window_specs)
    sm.plt.plots = Vector{StellarModels.JemsPlot}(undef, no_of_plots)
    for j = 1:no_of_plots
        this_axis = Axis(sm.plt.fig[sm.opt.plotting.window_layouts[j]...])
        this_type = Symbol(sm.opt.plotting.window_specs[j])
        sm.plt.plots[j] = StellarModels.JemsPlot(this_axis, this_type)
        if sm.plt.plots[j].type == :profile && length(sm.opt.plotting.profile_alt_yaxes) > 0
            sm.plt.plots[j].alt_ax = Axis(sm.plt.fig[sm.opt.plotting.window_layouts[j]...])
            sm.plt.plots[j].alt_ax.yaxisposition = :right
            linkxaxes!(sm.plt.plots[j].ax, sm.plt.plots[j].alt_ax)
        end
    end
end

"""
    create_HR_observables!(sm::StellarModel, plot::StellarModels.JemsPlot)

creates teff and L observables and adds them to the observable list of the given plot
"""
function create_HR_observables!(sm::StellarModel, plot::StellarModels.JemsPlot)
    plot.x_obs = Dict{Symbol, Observable}()
    teff = exp(sm.esi.lnT[sm.nz])
    plot.x_obs[:Teff_now] = Observable{Float64}(teff)
    plot.x_obs[:Teff] = Observable(Float64[])
    push!(plot.x_obs[:Teff][], teff)
    plot.y_obs = Dict{Symbol, Observable}()
    plot.y_obs[:L_now] = Observable{Float64}(sm.esi.L[sm.nz])
    plot.y_obs[:L] = Observable(Float64[])
    push!(plot.y_obs[:L][], sm.esi.L[sm.nz])
end

"""
    init_HR_plot!(ax::Axis, Teff::Observable, L::Observable, Teff_now::Observable, L_now::Observable;
                      line_kwargs=Dict(), scatter_kwargs=Dict())

Sets up the plot elements for an HRD
"""
function make_HR_plot!(ax::Axis, Teff::Observable, L::Observable, Teff_now::Observable, L_now::Observable;
                       line_kwargs=Dict(), scatter_kwargs=Dict())
    ax.xlabel = L"\log_{10}(T_\mathrm{eff}/[K])"
    ax.ylabel = L"L/L_\odot"
    ax.xreversed = true
    lines!(ax, Teff, L; line_kwargs...)
    scatter!(ax, Teff_now, L_now; scatter_kwargs...)
end

function create_profile_observables!(sm::StellarModel, plot::StellarModels.JemsPlot)
    xname = sm.opt.plotting.profile_xaxis
    xvals = StellarModels.profile_output_options[xname][2].((sm,), 1:(sm.nz))
    plot.x_obs = Dict{Symbol, Observable}()
    plot.x_obs[:profile_x] = Observable{Vector{Float64}}(xvals)

    ynames = sm.opt.plotting.profile_yaxes
    plot.y_obs = Dict{Symbol, Observable}()
    for name in ynames
        plot.y_obs[Symbol(name)] = 
            Observable{Vector{Float64}}(StellarModels.profile_output_options[name][2].((sm,), 1:(sm.nz)))
    end
    ynames = sm.opt.plotting.profile_alt_yaxes
    plot.alt_y_obs = Dict{Symbol,Observable}()
    for name in ynames
        plot.alt_y_obs[Symbol(name)] = Observable{Vector{Float64}}(StellarModels.profile_output_options[name][2].((sm,),
                                                                                                              1:(sm.nz)))
    end
end

"""
    function make_profile_plot!(ax::Axis, xvals::Observable, yvals::Observable, 
                            xlabel::AbstractString="", ylabels::Dict{Symbol, <:AbstractString}=Dict(); line_kwargs=Dict())

Plot a line for each entry in the 'yvals' observable, and put the given 'ylabels' in a legend.
"""
function make_profile_plot!(ax::Axis, xvals::Observable, yvals::Dict{Symbol,Observable},
                            xlabel::AbstractString="", ylabels::Dict{Symbol,<:AbstractString}=Dict();
                            alt_ax::Axis=nothing, alt_yvals::Dict{Symbol,Observable}=nothing, 
                            alt_ylabels::Dict{Symbol,<:AbstractString}=nothing,
                            line_kwargs=Dict())
    ax.xlabel = xlabel
    (color, state) = iterate(colors)
    for (name, yline) in yvals
        lines!(ax, xvals, yline, line_kwargs..., color=color, label=ylabels[name])
        (color, state) = iterate(colors, state)
    end
    axislegend(ax, position=:lt)
    if ! isnothing(alt_ax)
        for (name, yline) in alt_yvals
            lines!(alt_ax, xvals, yline, line_kwargs..., color=color, label=alt_ylabels[name])
        end
        axislegend(alt_ax)
    end

end

"""
    update_plots!(sm::StellarModel)

Updates all plots currently being displayed, by collecting appropriate data and notifying observables
"""
function update_plotting!(sm::StellarModel)
    for plot in sm.plt.plots
        if plot.type == :HR
            push!(plot.x_obs[:Teff][], exp(sm.esi.lnT[sm.nz]))
            plot.x_obs[:Teff_now].val = exp(sm.esi.lnT[sm.nz])
            push!(plot.y_obs[:L][], sm.esi.L[sm.nz])
            plot.y_obs[:L_now].val = sm.esi.L[sm.nz]
        elseif plot.type == :profile
            plot.x_obs[:profile_x].val = 
                StellarModels.profile_output_options[sm.opt.plotting.profile_xaxis][2].((sm,), 1:(sm.nz))
            for (key, obs) in pairs(plot.y_obs)
                obs.val = StellarModels.profile_output_options[String(key)][2].((sm,), 1:(sm.nz))
            end
            for (key, obs) in pairs(plot.alt_y_obs)
                obs.val = StellarModels.profile_output_options[String(key)][2].((sm,), 1:(sm.nz))
            end
        end
    end
    for plot in sm.plt.plots
        for xobs in values(plot.x_obs)
            notify(xobs)  # notifying only the x observables should replot everything
        end
        try
            autolimits!(plot.ax)
        catch  # catches Float32 errors in Makie
        end  # empty blocks allowed in julia!
        if !isnothing(plot.alt_ax)
            try
                autolimits!(plot.alt_ax)
            catch
            end
        end
    end
end

"""
    end_of_evolution(sm::StellarModel)

Perform end of evolution actions
"""
function end_of_evolution(sm::StellarModel)
    if sm.opt.plotting.wait_at_termination
        GLMakie.wait(sm.plt.scr)
    end
end

end  # end module Plotting
