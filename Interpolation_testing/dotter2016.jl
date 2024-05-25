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
gridpath = "DualRuns/DualGrid2"
path = "DualRuns/DualGrid/logM_-0.1_X_0.7381_.history.hdf5"
get_logM(path) = parse(Float64, split(split(path, "logM_")[2], "_")[1])

#get ALL filepaths in gridpath
historypaths = filter(x -> occursin(".history.hdf5", x), readdir(gridpath))
profilepaths = filter(x -> occursin(".profiles.hdf5", x), readdir(gridpath))
N = length(historypaths)
println("$N history files found")
models = Dict()
modeltracks = Dict()
X_dual         = ForwardDiff.Dual{}(0.7381,  0.0,1.0,0.0)
Z_dual         = ForwardDiff.Dual{}(0.0134,  0.0,0.0,1.0)
Dfraction_dual = ForwardDiff.Dual{}(0.000312,0.0,0.0,0.0)
R_dual         = ForwardDiff.Dual{}(100*RSUN,0.0,0.0,0.0)

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
    logM_dual = ForwardDiff.Dual{}(logM, 1.0,0.0,0.0)
    initial_params = [logM_dual, X_dual, Z_dual, Dfraction_dual, R_dual]
    model = de.Model_constructor(history_dual, profiles_dual, initial_params, inititial_params_names)
    track = nothing
    try 
        X_init = 0.4; X_end = 0.0000001
        X_init = 0.99*model.history_value.X_center[1]; X_end = 0.0001
        track = de.Track(model,X_init, X_end, 1000)
    catch 
        println(" Track FAILED for logM = $logM")
        continue
    end
    println(" Track OK for logM = $logM")
    modeltracks[logM] = track; models[logM] = model
    push!(all_logMs, logM)
end
##
fig = Figure();
ax = Axis(fig[1, 1],ylabel=L"\log (M/M_{\odot})", xlabel=L"\log(\text{Age/yr})")
ax.xgridvisible=true; ax.ygridvisible=true
for (logM, track) in modeltracks
    scatter!.(ax, log10.(models[logM].history_value.star_age),logM, color=:yellow,alpha=0.5,markersize=15)
    scatter!.(ax, log10.(track.history_value.star_age)       ,logM, color=:black)
end
ax2 = Axis(fig[1, 2],ylabel=L"\log (M/M_{\odot})", xlabel=L"Central hydrogen $X$",xreversed=true)
for (logM, track) in modeltracks
    scatter!.(ax2, models[logM].history_value.X_center,logM, color=:yellow,alpha=0.5,markersize=15)
    scatter!.(ax2, track.history_value.X_center       ,logM, color=:black)
end

fig
##
track = modeltracks[0.0];

##
fig = Figure();
#hrd
ax = Axis(fig[1, 1],ylabel=L"\log T_{\text{eff}}", xlabel=L"\log L/L_{\odot}",xreversed=true)
model = models[-0.45]
scatter!(ax, log10.(model.history_value.T_surf),log10.(model.history_value.L_surf), color=:yellow,alpha=0.5,markersize=15)
fig
##