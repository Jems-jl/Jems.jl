using DataFrames
using ForwardDiff, DataInterpolations
import ForwardDiff.Dual
using Jems.StellarModels
using Jems.DualExtrapolation; const de = DualExtrapolation
using Interpolations, Dierckx
using HDF5
using Jems.Constants
using DataInterpolations: CubicSpline
using CairoMakie, LaTeXStrings, MathTeXEngine, Makie.Colors, PlotUtils
basic_theme = Theme(fonts=(regular=texfont(:text), bold=texfont(:bold),
                           italic=texfont(:italic), bold_italic=texfont(:bolditalic)),
                    fontsize=30, size=(1000, 750), linewidth=7,
                    Axis=(xlabelsize=40, ylabelsize=40, titlesize=40, xgridvisible=false, ygridvisible=false,
                          spinewidth=2.5, xminorticksvisible=true, yminorticksvisible=true, xtickalign=1, ytickalign=1,
                          xminortickalign=1, yminortickalign=1, xticksize=14, xtickwidth=2.5, yticksize=14,
                          ytickwidth=2.5, xminorticksize=7, xminortickwidth=2.5, yminorticksize=7, yminortickwidth=2.5,
                          xticklabelsize=35, yticklabelsize=35, xticksmirrored=true, yticksmirrored=true),
                    Legend=(patchsize=(70, 10), framevisible=false, patchlabelgap=20, rowgap=10))
set_theme!(basic_theme)



##
gridpath = "DualRuns/DualGrid"
path = "DualRuns/DualGrid/logM_-0.1_X_0.7381_.history.hdf5"
get_logM(path) = parse(Float64, split(split(path, "logM_")[2], "_")[1])

#get ALL filepaths in gridpath
historypaths = filter(x -> occursin(".history.hdf5", x), readdir(gridpath))
profilepaths = filter(x -> occursin(".profiles.hdf5", x), readdir(gridpath))
N = length(historypaths)
println("$N history files found")
models = Dict()
modeltracks = Dict()
X_dual         = ForwardDiff.Dual{}(0.7381,  0.0,1.0,0.0,0.0,0.0)
Z_dual         = ForwardDiff.Dual{}(0.0134,  0.0,0.0,1.0,0.0,0.0)
Dfraction_dual = ForwardDiff.Dual{}(0.000312,0.0,0.0,0.0,1.0,0.0)
R_dual         = ForwardDiff.Dual{}(100*RSUN,0.0,0.0,0.0,0.0,1.0)

inititial_params_names = [:logM, :X, :Z, :Dfraction, :R]
all_logMs = []
i=8
only_use_these_logMs = [0.0,0.02]
for i in 1:N
    historypath = joinpath(gridpath, historypaths[i])
    profilepath = joinpath(gridpath, profilepaths[i])
    @assert get_logM(historypath) == get_logM(profilepath)
    logM = get_logM(historypath)
    if !isempty(only_use_these_logMs) if logM ∉ only_use_these_logMs; continue; end; end
    history_dual, profiles_dual = de.bookkeeping(historypath, profilepath,5)
    logM_dual = ForwardDiff.Dual{}(logM, 1.0,0.0,0.0,0.0,0.0)
    initial_params = [logM_dual, X_dual, Z_dual, Dfraction_dual, R_dual]
    model = de.Model_constructor(history_dual, profiles_dual, initial_params, inititial_params_names)
    models[logM] = model
    track = nothing
    try 
        X_init = 0.4; X_end = 0.0000001
        X_init = 0.99*model.history_value.X_center[1]; X_end = 0.001
        track = de.Track(model,X_init, X_end, 1000)
    catch 
        println(" Track FAILED for logM = $logM")
        continue
    end
    println(" Track OK for logM = $logM")
    modeltracks[logM] = track
    push!(all_logMs, logM)
end
## SMALL MASS CHANGES
extrapolGrid_00 = de.ExtrapolGrid(modeltracks[0.0], .- all_logMs);
extrapolGrid_00 = de.ExtrapolGrid(modeltracks[0.0], [0.01,0.02,-0.01,0.03]);
extrapolGrid_00 = de.ExtrapolGrid(modeltracks[0.0], [0.3]);
extrapolGrid_00 = de.ExtrapolGrid(modeltracks[0.0], [0.02]);
## 
fig = Figure(figsize=(2000,1000))
ax = Axis(fig[1,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$", title=L"Extrapolation from $\log M = 0.0$",xreversed=true)
de.plot!(modeltracks[0.0], ax ;scatter=false, label = L"JEMS $\log M = 0.0$ (interpolated)",color=:black)
scatter!(modeltracks[0.0].history_value.logT, modeltracks[0.0].history_value.logL;alpha=1.0, markersize=5,color=:white, label="JEMS Output")
#dM = 0.01; plot!(modeltracks[dM], ax; label = L"JEMS $\log M = %$dM $")
dM = 0.02; de.plot!(modeltracks[dM], ax; label = L"JEMS $\log M = %$dM $")
scatter!(modeltracks[0.02].history_value.logT, modeltracks[0.02].history_value.logL;alpha=1.0, markersize=5,color=:white, label="JEMS Output")
de.plot_arrows!(ax,modeltracks[0.0],0:0.01:1,dM)
de.plot!(extrapolGrid_00, ax; plot_original=false,scatter=true)
leg = Legend(fig[1,2],ax)
#xlims!(ax, [4.11,4.06]); ylims!(ax, [1.14,1.25])
fig 
##
keepx = ax.xaxis.attributes.limits.val; keepy = ax.yaxis.attributes.limits.val 
logM = 0.0 ;lines!(ax, models[logM].history_value.logT, models[logM].history_value.logL;alpha=0.1, color=:black, label="Original track")
logM = 0.02 ;lines!(ax, models[logM].history_value.logT, models[logM].history_value.logL;alpha=0.1, color=:black, label="Original track")
xlims!(ax, reverse(keepx)); ylims!(ax, keepy)
fig

## LARGER MASS CHANGES
extrapolGrid_00 = de.ExtrapolGrid(modeltracks[0.0], [0.1,0.2,0.3]);
fig = Figure(figsize=(2000,1500))
ax = Axis(fig[1,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$", title=L"Extrapolation from $\log M = 0.0$",xreversed=true)
de.plot!(modeltracks[0.0], ax ; label = L"JEMS $\log M = 0.0$",color=:black,scatter=false)
#de.plot!(modeltracks[0.05], ax; label = L"JEMS $\log M = 0.1$", color=:blue)
de.plot!(modeltracks[0.1], ax; label = L"JEMS $\log M = 0.1$", color=:blue)
de.plot!(modeltracks[0.2], ax; label = L"JEMS $\log M = 0.2$", color=:blue)
de.plot!(modeltracks[0.3], ax; label = L"JEMS $\log M = 0.3$", color=:blue)
de.plot!(extrapolGrid_00, ax; plot_original=false)
leg = Legend(fig[1,2],ax)
fig

## ZETA diagnostic PLOTS
extrapolGrid = de.ExtrapolGrid(modeltracks[0.0], all_logMs);
logMs_used = all_logMs
logMs_used = [0.01,0.1,0.3]
extrapolGrid = de.ExtrapolGrid(modeltracks[0.0],logMs_used);
##
fig = Figure(figsize=(3000,3000))
ax1 = Axis(fig[1,1], ylabel = L"$\log (L / L_\odot)$")
ax2 = Axis(fig[1,2], ylabel = L"$\log (T_{\text{eff}} / K)$")
scatter!(ax1, modeltracks[0.0].zeta, modeltracks[0.0].logL_val ; label = L"JEMS $\log M = 0.0$",color=:black)
scatter!(ax2, modeltracks[0.0].zeta, modeltracks[0.0].logT_val ; label = L"JEMS $\log M = 0.0$",color=:black)
for track in extrapolGrid.extrapoltracks
    logM = track.delta_logM + 0.0
    label = L"$\Delta \log M = %$logM $"
    lines!(ax1, track.zeta, track.logL_val, label = label, color=extrapolGrid.colors_dic[track],linewidth = 10)
    lines!(ax2, track.zeta, track.logT_val, label = label, color=extrapolGrid.colors_dic[track],linewidth = 10)
end
for logM in logMs_used #    actual JEMS tracks
    lines!(ax1, modeltracks[logM].zeta, modeltracks[logM].logL_val ,color=:black,linewidth=2)
    lines!(ax2, modeltracks[logM].zeta, modeltracks[logM].logT_val ,color=:black,linewidth=2)
end

ax3 = Axis(fig[2,1], xlabel = L"$\zeta$", ylabel = L"\partial \log L / \partial \log M")
hlines!(ax3, [3], color=:black, linestyle=:dot,alpha=0.4)
scatter!(ax3, modeltracks[0.0].zeta, modeltracks[0.0].logL_partial, label = "Original Track", color=:black)
ax4 = Axis(fig[2,2], xlabel = L"$\zeta$", ylabel = L"\partial \log T / \partial \log M")
scatter!(ax4, modeltracks[0.0].zeta, modeltracks[0.0].logT_partial, label = "Original Track", color=:black)

for logM in logMs_used
    lines!(ax3, modeltracks[logM].zeta, modeltracks[logM].logL_partial, label = "$logM",alpha=0.3,linewidth=5)
    lines!(ax4, modeltracks[logM].zeta, modeltracks[logM].logT_partial, label = "$logM",alpha=0.3,linewidth=5)
end

fig[1,3]   = Legend(fig, ax1)
fig[3,1:3] = Legend(fig, ax3,orientation=:horizontal,tellwidth=true)
fig
##

