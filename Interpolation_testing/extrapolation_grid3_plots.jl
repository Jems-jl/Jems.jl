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
gridpath = "DualRuns/DualGrid3"
get_logM(path) = parse(Float64, split(split(path, "logM_")[2], "_")[1])
#get ALL filepaths in gridpath
historypaths = filter(x -> occursin(".history.hdf5", x), readdir(gridpath))
profilepaths = filter(x -> occursin(".profiles.hdf5", x), readdir(gridpath))
N = length(historypaths)
println("$N history files found")
models = Dict()
modeltracks = Dict()
X_dual         = ForwardDiff.Dual{}(0.7154,  0.0,1.0,0.0)
Z_dual         = ForwardDiff.Dual{}(0.0142,  0.0,0.0,1.0)
Dfraction_dual = ForwardDiff.Dual{}(0.0,0.0,0.0,0.0)
R_dual         = ForwardDiff.Dual{}(100*RSUN,0.0,0.0,0.0)

inititial_params_names = [:logM, :X, :Z]
all_logMs = []
i=8
only_use_these_logMs = []#[0.0] #[0.0,0.02]
for i in 1:N
    historypath = joinpath(gridpath, historypaths[i])
    profilepath = joinpath(gridpath, profilepaths[i])
    @assert get_logM(historypath) == get_logM(profilepath)
    logM = get_logM(historypath)
    if !isempty(only_use_these_logMs) if logM ∉ only_use_these_logMs; continue; end; end
    history_dual, profiles_dual = de.bookkeeping(historypath, profilepath,3)
    logM_dual = ForwardDiff.Dual{}(logM, 1.0,0.0,0.0)
    initial_params = [logM_dual, X_dual, Z_dual]
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
end; println("Number of OK tracks: $(length(all_logMs))")
## SMALL MASS CHANGES
extrapolGrid_00 = de.ExtrapolGrid(modeltracks[0.0], [0.05,-0.05]);
##  EERST VOOR GRID2 GEDAAN, dan voor cubic en GRID3
fig = Figure(size=(1500,800))
ax = Axis(fig[1,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$", title=L"Extrapolation from $\log M = 0.0$",xreversed=true)
de.plot!(modeltracks[0.0], ax ;scatter=false, label = L"JEMS $\log M = 0.0$",color=:black,linewidth=10)
scatter!(modeltracks[0.0].history_value.logT, modeltracks[0.0].history_value.logL;alpha=1.0, markersize=10,color=:white)
#dM = 0.01; plot!(modeltracks[dM], ax; label = L"JEMS $\log M = %$dM $")
dM = 0.05; de.plot!(modeltracks[dM], ax; label = L"JEMS $\log M = %$dM $",scatter=false)
de.plot_arrows!(ax,modeltracks[0.0],0:0.01:1,dM; linewidth=0.3)

dM = -0.05; de.plot!(modeltracks[dM], ax; label = L"JEMS $\log M = %$dM $",color=:red,scatter=false)
de.plot_arrows!(ax,modeltracks[0.0],0:0.01:1,dM;linewidth=0.3)

de.plot!(extrapolGrid_00, ax; plot_original=false,scatter=true)
leg = Legend(fig[1,2],ax)
#xlims!(ax, [4.11,4.06]); ylims!(ax, [1.14,1.25])
fig 
##
figpath  = "Figures/extrapolation_firstExample_cubic.png"; save(figpath, fig; px_per_unit=3); @show figpath
##
keepx = ax.xaxis.attributes.limits.val; keepy = ax.yaxis.attributes.limits.val 
logM = 0.0 ;lines!(ax, models[logM].history_value.logT, models[logM].history_value.logL;alpha=0.1, color=:black, label="Original track")
logM = 0.02 ;lines!(ax, models[logM].history_value.logT, models[logM].history_value.logL;alpha=0.1, color=:black, label="Original track")
xlims!(ax, reverse(keepx)); ylims!(ax, keepy)
fig



extrapolGrid_00.extrapoltracks
extrapolGrid_00.delta_logM


## MASTERPLOT EXTRAPOLATION
logMs =  [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
extrapolGrid_00 = de.ExtrapolGrid(modeltracks[0.0],logMs);
fig = Figure(size=(800,600))
ax = Axis(fig[1,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L" $\log (L/L_\odot)$", title=L"Extrapolation from $\log M = 0.0$",xreversed=true)
ax.xlabelsize = 30; ax.ylabelsize = 30; ax.titlesize = 30
original = de.plot!(modeltracks[0.0], ax ; label = L"JEMS $\log M = 0.0$",color=:black,scatter=false,linewidth=4)
jemsline = nothing
for (delta,track) in zip(extrapolGrid_00.delta_logM, extrapolGrid_00.extrapoltracks)
    zams = (track.logT_val[1],track.logL_val[1])
    scatter!(ax, zams[1], zams[2], markersize=15, color=extrapolGrid_00.colors_dic[track])
    text!(ax, zams[1]+0.08, zams[2]-0.10, text = L"$ %$delta", color=:black, fontsize=20)
end

for logM in logMs
    try 
        jemsline = de.plot!(modeltracks[logM], ax; scatter=false,linewidth=2,color=:black) 
        scatter!(modeltracks[logM].history_value.logT[1], modeltracks[logM].history_value.logL[1];alpha=1.0, markersize=8,color=:black)
    catch
        println("Failed for logM = $logM")
    end
end
de.plot!(extrapolGrid_00, ax; scatter=false,linewidth=8,alpha=0.5)
#legend only containing info on jemsline
axislegend(ax, [jemsline, original], ["JEMS tracks", L"JEMS $\log M = 0.0$"])

de.plot_arrows!(ax,modeltracks[0.0],0:0.1:1,-1; linewidth=0.3)
de.plot_arrows!(ax,modeltracks[0.0],0:0.1:1,-0.8; linewidth=2,alpha=0.1)
fig

## EXTRA PLOT OMDRAAIING
logMs =  [-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1]#,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
extrapolGrid_00 = de.ExtrapolGrid(modeltracks[0.0],logMs);
fig = Figure(size=(1000,500))
ax = Axis(fig[1,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L" $\log (L/L_\odot)$",xreversed=true)
ax2 = Axis(fig[1,2], xlabel=L"$\log (T_{\text{eff}} / K)$", xreversed=true)
ax.xlabelsize = 30; ax.ylabelsize = 30; ax.titlesize = 30
original = de.plot!(modeltracks[0.0], ax ; label = L"JEMS $\log M = 0.0$",color=:black,scatter=false,linewidth=8)
jemsline = nothing
for (delta,track) in zip(extrapolGrid_00.delta_logM, extrapolGrid_00.extrapoltracks)
    zams = (track.logT_val[1],track.logL_val[1])
    scatter!(ax, zams[1], zams[2], markersize=15, color=extrapolGrid_00.colors_dic[track])
    scatter!(ax2, zams[1], zams[2], markersize=15, color=extrapolGrid_00.colors_dic[track])
    text!(ax, zams[1]+0.08, zams[2]-0.10, text = L"$ %$delta", color=:black, fontsize=20)
    text!(ax2, zams[1]+0.08, zams[2]-0.10, text = L"$ %$delta", color=:black, fontsize=20)
end

de.plot!(extrapolGrid_00, ax; plot_original=false, scatter=false,linewidth=8,alpha=0.5)
de.plot!(extrapolGrid_00, ax2; plot_original=false, scatter=false,linewidth=8,alpha=0.5)

for logM in logMs 
    de.plot_arrows!(ax,modeltracks[0.0],0.0:0.1:1,logM; linewidth=0.1) 
    de.plot_arrows!(ax2,modeltracks[0.0],0.0:0.01:1,logM; linewidth=0.2) 
end

de.plot_arrows!(ax, modeltracks[0.0],[0.0],-1; linewidth=5,color=:red)
de.plot_arrows!(ax2,modeltracks[0.0],[0.0],-1; linewidth=5,color=:red)

xlims!(ax2,3.77,3.67); ylims!(ax2,-2.05,-1.69)
text!(ax2, 3.749,-2.03, text = "-1.0")
ax2.xticks = [3.75,3.70]
fig
##
figpath  = "Figures/extrapolation00_masterplot_omdraaiing.png"; save(figpath, fig; px_per_unit=3); @show figpath
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

