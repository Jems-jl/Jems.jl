module Plotting

using GLMakie, LaTeXStrings, MathTeXEngine, Jems.StellarModels

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

    init_figure!(sm)
    sm.plt.obs = Dict()
    sm.plt.axno = Dict()

    n = 1
    if "HR" in sm.opt.plotting.window_flags
        create_HR_observables!(sm)
        init_HR_plot!(sm.plt.axs[n], sm.plt.obs[:Teff], sm.plt.obs[:L], sm.plt.obs[:Teff_now],
                    sm.plt.obs[:L_now]; scatter_kwargs=Dict(:color => "red"))
        sm.plt.axno[:HR] = n
        n += 1
    end
    if "profile" in sm.opt.plotting.window_flags
        create_profile_observables!(sm)
        init_profile_plot!(sm.plt.axs[n], sm.plt.obs[:profile_x], sm.plt.obs[:profile_y], "xlabel", "ylabel")
        sm.plt.axno[:profile] = n
        n += 1
    end
    sm.plt.scr = display(sm.plt.fig)
end

"""
    init_figure!(sm::StellarModel)

Initializes the plotting figure and adds axes
"""
function init_figure!(sm::StellarModel)
    sm.plt.fig = Figure()
    no_of_plots = determine_no_plots(sm)
    sm.plt.axs = Vector{Makie.Block}(undef, no_of_plots)
    for i=1:no_of_plots
        sm.plt.axs[i] = Axis(sm.plt.fig[i, 1])
    end
end

function determine_no_plots(sm::StellarModel)
    return length(sm.opt.plotting.window_flags)
end

"""
    create_HR_observables!(sm::StellarModel)

creates teff and L observables and adds them to the observable list
"""
function create_HR_observables!(sm::StellarModel)
    teff = exp(sm.esi.lnT[sm.nz])
    sm.plt.obs[:Teff_now] = Observable{Float64}(teff)
    sm.plt.obs[:Teff] = Observable(Float64[])
    push!(sm.plt.obs[:Teff][], teff)
    sm.plt.obs[:L_now] = Observable{Float64}(sm.esi.L[sm.nz])
    sm.plt.obs[:L] = Observable(Float64[])
    push!(sm.plt.obs[:L][], sm.esi.L[sm.nz])
end

"""

    init_HR_plot!(ax::Axis, Teff::Observable, L::Observable, Teff_now::Observable, L_now::Observable;
                      line_kwargs=Dict(), scatter_kwargs=Dict())

Sets up the plot elements for an HRD
"""
function init_HR_plot!(ax::Axis, Teff::Observable, L::Observable, Teff_now::Observable, L_now::Observable;
                      line_kwargs=Dict(), scatter_kwargs=Dict())
    ax.xlabel = L"\log_{10}(T_\mathrm{eff}/[K])"
    ax.ylabel = L"\log_{10}(L/L_\odot)"
    ax.xreversed = true
    lines!(ax, Teff, L; line_kwargs...)
    scatter!(ax, Teff_now, L_now; scatter_kwargs...)
end


function create_profile_observables!(sm::StellarModel)
    xname = sm.opt.plotting.profile_xaxis
    yname = sm.opt.plotting.profile_yaxis
    xvals = StellarModels.profile_output_options[xname][2].((sm,), 1:sm.nz)
    yvals = StellarModels.profile_output_options[yname][2].((sm,), 1:sm.nz)
    sm.plt.obs[:profile_x] = Observable{Vector{Float64}}(xvals)
    sm.plt.obs[:profile_y] = Observable{Vector{Float64}}(yvals)
end

function init_profile_plot!(ax::Axis, xvals::Observable, yvals::Observable, xlabel::String="", ylabel::String="", line_kwargs=Dict())
    ax.xlabel = xlabel
    ax.ylabel = ylabel
    lines!(ax, xvals, yvals, line_kwargs...)
end

"""
    update_plots!(sm::StellarModel)

Updates all plots currently being displayed, by collecting appropriate data and notifying observables
"""
function update_plots!(sm::StellarModel)
    if "HR" in sm.opt.plotting.window_flags
        push!(sm.plt.obs[:Teff][], exp(sm.esi.lnT[sm.nz]))
        sm.plt.obs[:Teff_now][] = exp(sm.esi.lnT[sm.nz])
        push!(sm.plt.obs[:L][], sm.esi.L[sm.nz])
        sm.plt.obs[:L_now][] = sm.esi.L[sm.nz]
        notify(sm.plt.obs[:Teff])  # needed because we change these vars in place
        autolimits!(sm.plt.axs[sm.plt.axno[:HR]])
    end
    if "profile" in sm.opt.plotting.window_flags
        sm.plt.obs[:profile_x].val = StellarModels.profile_output_options[sm.opt.plotting.profile_xaxis][2].((sm,), 1:sm.nz)
        sm.plt.obs[:profile_y][] = StellarModels.profile_output_options[sm.opt.plotting.profile_yaxis][2].((sm,), 1:sm.nz)
    end
end

"""
    end_of_evolution(sm::StellarModel)

Perform end of evolution actions
"""
function end_of_evolution(sm::StellarModel)
    GLMakie.wait(sm.plt.scr)
end

end
