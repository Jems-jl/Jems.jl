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
gridpath = "DualRuns/DualGrid2"
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
only_use_these_logMs = []
for i in 1:N
    historypath = joinpath(gridpath, historypaths[i])
    profilepath = joinpath(gridpath, profilepaths[i])
    @assert get_logM(historypath) == get_logM(profilepath)
    logM = get_logM(historypath)
    if !isempty(only_use_these_logMs) if logM ∉ only_use_these_logMs; continue; end; end
    history_dual, profiles_dual = de.bookkeeping(historypath, profilepath,3)
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
##
interpol_asked = [0.05,0.1,0.15]
interpolTrackxs = Dict()
modeltrack1, modeltrack2 = modeltracks[0.0], modeltracks[0.2]
extrapolGrid1 = de.ExtrapolGrid(modeltrack1, interpol_asked .- modeltrack1.logM);
extrapolGrid2 = de.ExtrapolGrid(modeltrack2, interpol_asked .- modeltrack2.logM);
for logM in interpol_asked
    interpolTrack = de.InterpolTrack(modeltrack1, modeltrack2,logM);
    interpolTrackxs[logM] = interpolTrack
end
##
fig = Figure(size=(2000,1000))
ax = Axis(fig[1,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$", title=L"Interpolation between $\log M = 0$ and $0.2",xreversed=true)
de.plot!(modeltrack1, ax; color=:red, scatter=false, label = "JEMS log M = " * string(modeltrack1.logM),linestyle=:dash, linewidth=5)
de.plot!(modeltrack2, ax; color=:blue, scatter=false,label = "JEMS log M = " * string(modeltrack2.logM),linestyle=:dot, linewidth=5)
for logM in interpol_asked
    interpolTrack = interpolTrackxs[logM]
    de.plot!(interpolTrack, ax; color=:lightblue, scatter=true,linewidth=5)
    de.plot!(modeltracks[logM], ax; color=:black, scatter=false,linewidth=2)
end
for extrapolTrack in extrapolGrid1.extrapoltracks
    de.plot!(extrapolTrack, ax; color=:red, scatter=false,linewidth=2,linestyle = :dash)
end
for extrapolTrack in extrapolGrid2.extrapoltracks
    de.plot!(extrapolTrack, ax; color=:blue, scatter=false,linewidth=2, linestyle = :dot)
end
axislegend(ax,position = :lb)
ax2 = Axis(fig[1,2], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$", title=L"Interpolation of $\log M = 0.05,0.10,0.15$",xreversed=true)
for logM in interpol_asked
    interpolTrack = interpolTrackxs[logM]
    de.plot!(interpolTrack, ax2; color=:lightblue, scatter=true,label="Interpol",linewidth=5)
    de.plot!(modeltracks[logM], ax2; color=:black, scatter=false,label="JEMS",linewidth=2)
end
for extrapolTrack in extrapolGrid1.extrapoltracks
    de.plot!(extrapolTrack, ax2; color=:red, scatter=false,label="Extrapol from below",linewidth=2,linestyle = :dash)
end
for extrapolTrack in extrapolGrid2.extrapoltracks
    de.plot!(extrapolTrack, ax2; color=:blue, scatter=false,label="Extrapol from above",linewidth=2, linestyle = :dot)
end
leg = Legend(fig[1,3],ax2;tellwidth=true)
fig
##



##
interpol_asked = [0.0]
interpolTrackxs = Dict()
modeltrack1, modeltrack2 = modeltracks[-0.1], modeltracks[0.1]
extrapolGrid1 = de.ExtrapolGrid(modeltrack1, interpol_asked .- modeltrack1.logM);
extrapolGrid2 = de.ExtrapolGrid(modeltrack2, interpol_asked .- modeltrack2.logM);
for logM in interpol_asked
    interpolTrack = de.InterpolTrack(modeltrack1, modeltrack2,logM);
    interpolTrackxs[logM] = interpolTrack
end
##
fig = Figure(size=(2000,1000))
ax = Axis(fig[1,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$", title=L"Interpolation between $\log M = -0.1$ and $0.1",xreversed=true)
de.plot!(modeltrack1, ax; color=:red, scatter=false, label = "JEMS log M = " * string(modeltrack1.logM),linestyle=:dash, linewidth=5)
de.plot!(modeltrack2, ax; color=:blue, scatter=false,label = "JEMS log M = " * string(modeltrack2.logM),linestyle=:dot, linewidth=5)
for logM in interpol_asked
    interpolTrack = interpolTrackxs[logM]
    de.plot!(interpolTrack, ax; color=:lightblue, scatter=true,linewidth=5)
    de.plot!(modeltracks[logM], ax; color=:black, scatter=false,linewidth=2)
end
for extrapolTrack in extrapolGrid1.extrapoltracks
    de.plot!(extrapolTrack, ax; color=:red, scatter=false,linewidth=2,linestyle = :dash)
end
for extrapolTrack in extrapolGrid2.extrapoltracks
    de.plot!(extrapolTrack, ax; color=:blue, scatter=false,linewidth=2, linestyle = :dot)
end
de.plot_arrows!(ax, modeltrack1, 0:0.1:1,0.1)
de.plot_arrows!(ax, modeltrack2, 0:0.1:1,-0.1)
axislegend(ax,position = :lb)
ax2 = Axis(fig[1,2], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$", title=L"Interpolation of $\log M =0.0$",xreversed=true)
for logM in interpol_asked
    interpolTrack = interpolTrackxs[logM]
    de.plot!(interpolTrack, ax2; color=:lightblue, scatter=true,label="Interpol",linewidth=5)
    #de.plot!(modeltracks[logM], ax2; color=:black, scatter=false,label="JEMS",linewidth=2)
end
for extrapolTrack in extrapolGrid1.extrapoltracks
    de.plot!(extrapolTrack, ax2; color=:red, scatter=false,label="Extrapol from below",linewidth=2,linestyle = :dash)
end
for extrapolTrack in extrapolGrid2.extrapoltracks
    de.plot!(extrapolTrack, ax2; color=:blue, scatter=false,label="Extrapol from above",linewidth=2, linestyle = :dot)
end
leg = Legend(fig[1,3],ax2;tellwidth=true)
fig


##
interpol_asked = [0.01]
interpolTrackxs = Dict()
modeltrack1, modeltrack2 = modeltracks[0.0], modeltracks[0.02]
extrapolGrid1 = de.ExtrapolGrid(modeltrack1, interpol_asked .- modeltrack1.logM);
extrapolGrid2 = de.ExtrapolGrid(modeltrack2, interpol_asked .- modeltrack2.logM);
for logM in interpol_asked
    interpolTrack = de.InterpolTrack(modeltrack1, modeltrack2,logM);
    interpolTrackxs[logM] = interpolTrack
end
##
fig = Figure(size=(2000,1000))
ax = Axis(fig[1,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$", title=L"Interpolation between $\log M = 0.0$ and $0.02",xreversed=true)
de.plot!(modeltrack1, ax; color=:red, scatter=false, label = "JEMS log M = " * string(modeltrack1.logM),linestyle=:dash, linewidth=5)
de.plot!(modeltrack2, ax; color=:blue, scatter=false,label = "JEMS log M = " * string(modeltrack2.logM),linestyle=:dot, linewidth=5)
for logM in interpol_asked
    interpolTrack = interpolTrackxs[logM]
    de.plot!(interpolTrack, ax; color=:lightblue, scatter=true,linewidth=5)
    de.plot!(modeltracks[logM], ax; color=:black, scatter=false,linewidth=2)
end
for extrapolTrack in extrapolGrid1.extrapoltracks
    de.plot!(extrapolTrack, ax; color=:red, scatter=false,linewidth=2,linestyle = :dash)
end
for extrapolTrack in extrapolGrid2.extrapoltracks
    de.plot!(extrapolTrack, ax; color=:blue, scatter=false,linewidth=2, linestyle = :dot)
end
de.plot_arrows!(ax, modeltrack1, 0:0.01:1,0.01)
de.plot_arrows!(ax, modeltrack2, 0:0.01:1,-0.01)
axislegend(ax,position = :lb)
ax2 = Axis(fig[1,2], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$", title=L"Interpolation of $\log M =0.01$",xreversed=true)
for logM in interpol_asked
    interpolTrack = interpolTrackxs[logM]
    de.plot!(interpolTrack, ax2; color=:lightblue, scatter=true,label="Interpol",linewidth=5)
    de.plot!(modeltracks[logM], ax2; color=:black, scatter=false,label="JEMS",linewidth=1)
end
for extrapolTrack in extrapolGrid1.extrapoltracks
    de.plot!(extrapolTrack, ax2; color=:red, scatter=false,label="Extrapol from below",linewidth=2,linestyle = :dash)
end
for extrapolTrack in extrapolGrid2.extrapoltracks
    de.plot!(extrapolTrack, ax2; color=:blue, scatter=false,label="Extrapol from above",linewidth=2, linestyle = :dot)
end
leg = Legend(fig[1,3],ax2;tellwidth=true)
fig