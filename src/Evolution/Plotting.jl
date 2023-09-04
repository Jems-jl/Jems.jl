
using GLMakie

function init_figure(sm::StellarModel)
    sm.plt.fig = Figure()
    sm.plt.axs = Vector{Makie.Block}(undef, 1)
    sm.plt.axs[1] = Axis(sm.plt.fig[1, 1])
end

function init_plot(fig::Figure, ax::Axis, o1::Observable, o2::Observable)
    scatter!(ax, o1, o2)
    return display(fig)
end

function set_profile_observables(ax::Axis, profile, x::String, y::String)
    X = @lift($profile[!, x])
    Y = @lift($profile[!, y])
    model_number_str = @lift("model number=$(parse(Int,$pname))")

    profile_line = lines!(ax, X, Y; label="real profile")
    profile_text = text!(ax, 0.7, 0.0; text=model_number_str)
end

function set_history_observables(ax::Axis, sm::StellarModel)
    ax.xlabel = L"\log_{10}(T_\mathrm{eff}/[K])"
    ax.ylabel = L"\log_{10}(L/L_\odot)"
    ax.xreversed = true
    xx = exp(sm.esi.lnT[sm.nz])
    sm.plt.obs = Dict()
    sm.plt.obs[:Teff] = Observable(Float64[])
    push!(sm.plt.obs[:Teff][], xx)
    sm.plt.obs[:L] = Observable(Float64[])
    push!(sm.plt.obs[:L][], sm.esi.L[sm.nz])
end

function push_observable!(o::Observable{Vector{TT}}, x::TT) where {TT<:Real}
    push!(o[], x)
end
