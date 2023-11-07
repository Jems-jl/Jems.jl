module Plotting

using GLMakie, LaTeXStrings, MathTeXEngine, Jems.StellarModels

basic_theme = Theme(fonts=(regular=texfont(:text), bold=texfont(:bold),
                           italic=texfont(:italic), bold_italic=texfont(:bolditalic)),
                    fontsize=30, resolution=(1000, 750), linewidth=7,
                    Axis=(xlabelsize=40, ylabelsize=40, titlesize=40, xgridvisible=false, ygridvisible=false,
                          spinewidth=2.5, xminorticksvisible=true, yminorticksvisible=true, xtickalign=1, ytickalign=1,
                          xminortickalign=1, yminortickalign=1, xticksize=14, xtickwidth=2.5, yticksize=14,
                          ytickwidth=2.5, xminorticksize=7, xminortickwidth=2.5, yminorticksize=7, yminortickwidth=2.5,
                          xticklabelsize=35, yticklabelsize=35, xticksmirrored=true, yticksmirrored=true),
                    Legend=(patchsize=(70, 10), framevisible=false, patchlabelgap=20, rowgap=10))
set_theme!(basic_theme)

GLMakie.set_window_config!(; float=true)  # place windows on top

"""
    init_plots(sm::StellarModel)

Sets up all observables to be traced this run, and creates the figure and axis where they will be plotted
"""
function init_plots!(sm::StellarModel)
    init_figure!(sm)
    sm.plt.obs = Dict()
    set_HR_observables!(sm)
    init_HR_plot!(sm.plt.axs[1], sm.plt.obs[:Teff], sm.plt.obs[:L], sm.plt.obs[:Teff_now],
                 sm.plt.obs[:L_now]; scatter_kwargs=Dict(:color => "red"))
    sm.plt.scr = display(sm.plt.fig)
end

function init_figure!(sm::StellarModel)
    sm.plt.fig = Figure()
    sm.plt.axs = Vector{Makie.Block}(undef, 1)
    sm.plt.axs[1] = Axis(sm.plt.fig[1, 1])
end

function init_HR_plot!(ax::Axis, Teff::Observable, L::Observable, Teff_now::Observable, L_now::Observable;
                      line_kwargs=Dict(), scatter_kwargs=Dict())
    ax.xlabel = L"\log_{10}(T_\mathrm{eff}/[K])"
    ax.ylabel = L"\log_{10}(L/L_\odot)"
    ax.xreversed = true
    lines!(ax, Teff, L; line_kwargs...)
    scatter!(ax, Teff_now, L_now; scatter_kwargs...)
end

function set_HR_observables!(sm::StellarModel)
    teff = exp(sm.esi.lnT[sm.nz])
    sm.plt.obs[:Teff_now] = Observable{Float64}(teff)
    sm.plt.obs[:Teff] = Observable(Float64[])
    push!(sm.plt.obs[:Teff][], teff)
    sm.plt.obs[:L_now] = Observable{Float64}(sm.esi.L[sm.nz])
    sm.plt.obs[:L] = Observable(Float64[])
    push!(sm.plt.obs[:L][], sm.esi.L[sm.nz])
end

function update_plots!(sm::StellarModel)
    push!(sm.plt.obs[:Teff][], exp(sm.esi.lnT[sm.nz]))
    sm.plt.obs[:Teff_now][] = exp(sm.esi.lnT[sm.nz])
    push!(sm.plt.obs[:L][], sm.esi.L[sm.nz])
    sm.plt.obs[:L_now][] = sm.esi.L[sm.nz]
    reset_limits!(sm.plt.axs[1])
    notify_observables(sm.plt.obs)
end

function notify_observables(os::Dict{Symbol,Observable})
    for o in values(os)
        notify(o)
    end
end

end
