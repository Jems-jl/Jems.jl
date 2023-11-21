"""
    create_HR_observables!(sm::StellarModel, plot::StellarModels.JemsPlot)

Creates teff and L observables and adds them to the observable list of the given plot
"""
function create_HR_observables!(plot::StellarModels.JemsPlot, sm::StellarModel)
    teff = exp(sm.esi.lnT[sm.nz])
    plot.x_obs[:Teff_now] = Observable{Float64}(teff)
    plot.x_obs[:Teff] = Observable(Float64[])
    push!(plot.x_obs[:Teff][], teff)
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

"""
    function update_HR_plot!(plot::StellarModels.JemsPlot, sm::StellarModel)

Updates the given `plot` with the relevant HR data from the stellar model `sm`.
"""
function update_HR_plot!(plot::StellarModels.JemsPlot, sm::StellarModel)
    push!(plot.x_obs[:Teff].val, exp(sm.esi.lnT[sm.nz]))
    plot.x_obs[:Teff_now].val = exp(sm.esi.lnT[sm.nz])
    push!(plot.y_obs[:L].val, sm.esi.L[sm.nz])
    plot.y_obs[:L_now].val = sm.esi.L[sm.nz]
end
