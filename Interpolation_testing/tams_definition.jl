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
function construct(X_tams)
models = Dict()
modeltracks = Dict()
X_dual         = ForwardDiff.Dual{}(0.7154,  0.0,1.0,0.0)
Z_dual         = ForwardDiff.Dual{}(0.0142,  0.0,0.0,1.0)

inititial_params_names = [:logM, :X, :Z]
all_logMs = []
i=8
only_use_these_logMs = [] #[0.0,0.02]
only_use_these_logMs = [0.0,0.05,-0.05] #[0.0,0.02]

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
        X_init = 0.99*model.history_value.X_center[1]; X_end = X_tams
        track = de.Track(model,X_init, X_end, 1000)
    catch 
        println(" Track FAILED for logM = $logM")
        continue
    end
    println(" Track OK for logM = $logM")
    modeltracks[logM] = track
    push!(all_logMs, logM)
end; println("Number of OK tracks: $(length(all_logMs))")
return models, modeltracks, all_logMs
end

X_tams_array = [0.01,0.001,0.0005,0.0001]
models_array = []
modeltracks_array = []
all_logMs_array = []

for X_tams in X_tams_array
    println("X_tams = $X_tams")
    models, modeltracks, all_logMs = construct(X_tams)
    push!(models_array, models)
    push!(modeltracks_array, modeltracks)
    push!(all_logMs_array, all_logMs)
end

## extrapolgrids for all xtams
deltas = [0.05,-0.05]
extrapolgrids = []
for modeltrack in modeltracks_array
    push!(extrapolgrids,de.ExtrapolGrid(modeltrack[0.0], deltas))
end;
## 
fig = Figure(size=(1500,1500))
ax1 = Axis(fig[1,1],                                        ylabel=L"$\log (L/L_\odot)$",  xreversed=true)
ax2 = Axis(fig[1,2],                                                                       xreversed=true)
ax3 = Axis(fig[2,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$",  xreversed=true)
ax4 = Axis(fig[2,2], xlabel=L"$\log (T_{\text{eff}} / K)$",                                xreversed=true)
axes = [ax1,ax2,ax3,ax4]
ax2.yticklabelsvisible = false; ax4.yticklabelsvisible = false
ax1.xticklabelsvisible = false; ax2.xticklabelsvisible = false
for (xtams, modeltrack, extrapolgrid, ax) in zip(X_tams_array,modeltracks_array, extrapolgrids, axes)
    de.plot!(modeltrack[0.0], ax ;scatter=false, label = L"JEMS $\log M = 0.0$",color=:black,linewidth=10)
    scatter!(ax,modeltrack[0.0].history_value.logT, modeltrack[0.0].history_value.logL;alpha=1.0, markersize=10,color=:white)
    text!(ax, 4.080, 1.37; text = L"$X_\text{TAMS} = %$xtams", color=:black, fontsize=40)
    dM = 0.05; de.plot!(modeltrack[dM], ax; label = L"JEMS $\log M = %$dM $",scatter=false)
    de.plot_arrows!(ax,modeltrack[0.0],0.9:0.01:1,dM; linewidth=0.3)

    dM = -0.05; de.plot!(modeltrack[dM], ax; label = L"JEMS $\log M = %$dM $",color=:red,scatter=false)
    de.plot_arrows!(ax,modeltrack[0.0],0.9:0.01:1,dM;linewidth=0.3)

    de.plot!(extrapolgrid, ax; plot_original=false,scatter=true)
end
linkxaxes!(axes...); linkyaxes!(axes...)
fig 
##
for (modeltracks, models, ax) in zip(modeltracks_array, models_array,axes)
    keepx = ax.xaxis.attributes.limits.val; keepy = ax.yaxis.attributes.limits.val 
    logM = 0.0 ;lines!(ax, models[logM].history_value.logT, models[logM].history_value.logL;alpha=0.05, color=:black, label="Original track")
    xlims!(ax, reverse(keepx)); ylims!(ax, keepy)
end; fig
##
figpath  = "Figures/extrapol_tams_definition.png"; save(figpath, fig; px_per_unit=3); @show figpath
##
fig = Figure(size=(1000,600))
axA = Axis(fig[1,1], ylabel=L"$\log L/L_\odot")
axB = Axis(fig[1,2], ylabel=L"$\log T/K")
axi = Axis(fig[2,1],  ylabel=L"$\frac{\partial \log L}{\partial \zeta}")
axii =Axis(fig[2,2], ylabel=L"$\frac{\partial \log T}{\partial \zeta}")
ax1 = Axis(fig[3,1],xlabel=L"\zeta", ylabel=L"$\frac{\partial \log L}{\partial \log M}")
ax2 = Axis(fig[3,2],xlabel=L"\zeta", ylabel=L"$\frac{\partial \log T}{\partial \log M}")
axA.xticklabelsvisible = false; axB.xticklabelsvisible = false; axi.xticklabelsvisible = false; axii.xticklabelsvisible = false
axA.ylabelsize = 30; axB.ylabelsize = 30; axi.ylabelsize = 30; axii.ylabelsize = 30; ax1.ylabelsize = 30; ax2.ylabelsize = 30
linkxaxes!(axA,ax1,axi,axB,ax2,axii)
colors = [:olivedrab1, :lime, :olive, :darkgreen]
colors = [:lightskyblue, :slateblue1, :skyblue4, :navyblue]
for (xtams, modeltracks, models,color) in zip(X_tams_array,modeltracks_array, models_array,colors)
    lines!(axA, modeltracks[0.0].zeta, modeltracks[0.0].logL_val ; label = L"$ %$xtams $",color=color,colormap=:viridis,linewidth=10)
    scatter!(axA, modeltracks[0.0].history_value.zeta, modeltracks[0.0].history_value.logL;alpha=1.0, markersize=10,color=:white)
    der = derivative_numerical(modeltracks[0.0].zeta, modeltracks[0.0].logL_val)
    scatter!(axi, modeltracks[0.0].zeta, der ; label = L"$ %$xtams $",color=color,colormap=:viridis)
end
for (xtams, modeltracks, models,color) in zip(X_tams_array,modeltracks_array, models_array,colors)
    lines!(axB, modeltracks[0.0].zeta, modeltracks[0.0].logT_val ; label = L"$ %$xtams $",color=color,colormap=:viridis,linewidth=10)
    scatter!(axB, modeltracks[0.0].history_value.zeta, modeltracks[0.0].history_value.logT;alpha=1.0, markersize=10,color=:white)
    der = derivative_numerical(modeltracks[0.0].zeta, modeltracks[0.0].logT_val)
    scatter!(axii, modeltracks[0.0].zeta, der ; label = L"$ %$xtams $",color=color,colormap=:viridis)
end
for (xtams, modeltracks, models,color) in zip(X_tams_array,modeltracks_array, models_array,colors)
    scatter!(ax1, modeltracks[0.0].zeta, modeltracks[0.0].logL_partial ; label = L"$ %$xtams $",color=color,colormap=:viridis)
end
for (xtams, modeltracks, models,color) in zip(X_tams_array,modeltracks_array, models_array,colors)
    scatter!(ax2, modeltracks[0.0].zeta, modeltracks[0.0].logT_partial ; label = L"$ %$xtams $",color=color,colormap=:viridis)
end
leg = Legend(fig[0,1:2],axA,orientation=:horizontal)
xlims!(ax1,0.9,1.1 ); xlims!(ax2,0.7,1.02)
ylims!(axA, 1.15,1.22); ylims!(axB, 4.06,4.1)
axA.yticksmirrored = false; axB.yticksmirrored = false; axi.yticksmirrored = false; axii.yticksmirrored = false; ax1.yticksmirrored = false; ax2.yticksmirrored = false; 
fig
##
savepath = "Figures/extrapol_tams_definition_zeta_diagnostic.png"; save(savepath, fig; px_per_unit=3); @show savepath
##

function derivative_numerical(x, y)
    N = length(x)
    dx = x[2] - x[1]
    dydx = zeros(N)
    dydx[1] = (y[2] - y[1]) / dx
    dydx[N] = (y[N] - y[N-1]) / dx
    for i in 2:N-1
        dydx[i] = (y[i+1] - y[i-1]) / (2dx)
    end
    return dydx
end

f(x) = x^2 + 5
f(Dual(5,1))