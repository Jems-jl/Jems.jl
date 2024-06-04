using DataFrames
using ForwardDiff, DataInterpolations
import ForwardDiff.Dual
using Jems.StellarModels
using Jems.DualExtrapolation; const de = DualExtrapolation
using Interpolations, Dierckx
using HDF5
using Jems.Constants, Statistics
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
struct Difference_with_JEMS
    logM
    diff_logL
    diff_logT
    diff_pyth
    diff_tot
    zeta
end
function Difference_with_JEMS(interpolTrack::de.InterpolTrack, modeltrack::de.Track)
    if interpolTrack.logM != modeltrack.logM
       @warn "logM of modeltrack (=$modeltrack.logM) and interpoltrack (=$interpolTrack.logM) not the same"
    end
    if length(interpolTrack.zeta) != length(modeltrack.zeta)
        throw(ArgumentError("Not the same zeta sampling"))
    end
    diff_logL = modeltrack.logL_val - interpolTrack.logL_val; logL_norm = maximum(abs.(diff_logL)) 
    diff_logT = modeltrack.logT_val - interpolTrack.logT_val; logT_norm = maximum(abs.(diff_logT))
    diff_pyth = sqrt.((diff_logL/logL_norm).^2 .+ (diff_logT/logT_norm).^2); 
    diff_tot = sum(diff_pyth)
    weight_L = 1 / abs(modeltrack.logL_val[end] - modeltrack.logL_val[1])
    weight_T = 1 / abs(modeltrack.logT_val[end] - modeltrack.logT_val[1])
    factor_L = weight_L / (weight_L + weight_T); factor_T = weight_T / (weight_L + weight_T)
    diff_pyth = sqrt.((diff_logL*factor_L).^2 .+ (diff_logT*factor_T).^2); diff_tot = sum(diff_pyth)
    return Difference_with_JEMS(interpolTrack.logM, diff_logL, diff_logT, diff_pyth, diff_tot, interpolTrack.zeta)
end
##
gridpath = "DualRuns/DualGrid"
gridpath = "DualRuns/DualGrid2"
gridpath = "DualRuns/DualGrid3"

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
        X_init = 0.99*model.history_value.X_center[1]; X_end = 0.01
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
interpol_asked = [0.02,0.03,0.04]
interpol_asked = [0.05,0.1,0.15]
interpol_asked = 0.05:0.05:0.55
interpol_asked = -0.95:0.05:-0.5
interpol_asked = -0.95:0.05:-0.65
interpolTrackxs = Dict{Float64,de.InterpolTrack}(); interpolTracks_dotter = Dict{Float64, de.InterpolTrack }()
modeltrack1, modeltrack2 = modeltracks[-0.1], modeltracks[0.6]
modeltrack1, modeltrack2 = modeltracks[-1.0], modeltracks[-0.6]
#extrapolGrid1 = de.ExtrapolGrid(modeltrack1, interpol_asked .- modeltrack1.logM);
#extrapolGrid2 = de.ExtrapolGrid(modeltrack2, interpol_asked .- modeltrack2.logM);

for logM in interpol_asked
    interpolTrack = de.InterpolTrack(modeltrack1, modeltrack2,logM);
    interpolTrack_dotter = de.InterpolTrack_dotter(modeltrack1, modeltrack2,logM);
    interpolTrackxs[logM] = interpolTrack; interpolTracks_dotter[logM] = interpolTrack_dotter
end
## DOTTER COMPARISON THESIS
fig = Figure(size=(1200,900))
ax = Axis(fig[1:2,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$",xreversed=true)
down = de.plot!(modeltrack1, ax; color=:red, scatter=false,linestyle=:dash, linewidth=5)
up = de.plot!(modeltrack2, ax; color=:blue, scatter=false,linestyle=:dot, linewidth=5)
jems_compare = nothing; interpol = nothing; dotter = nothing
for logM in interpol_asked
    interpolTrack = interpolTrackxs[logM]
    interpol = de.plot!(interpolTrack, ax; color=:lightgreen, scatter=false,linewidth=5)
    jems_compare = de.plot!(modeltracks[logM], ax; color=:black, scatter=false,linewidth=2)
    text!(ax, interpolTrack.logT_val[1]+0.04,interpolTrack.logL_val[1]-0.03; text = string(logM),fontsize=40)
    interpolTrack_dotter = interpolTracks_dotter[logM]
    dotter = de.plot!(interpolTrack_dotter, ax; color=:purple, scatter=false,linewidth=5,linestyle=:dashdot)
end

axislegend(ax,[up, jems_compare, interpol, dotter, down],["JEMS log M = " * string(modeltrack2.logM),"JEMS tracks for comparison", "Interpolated tracks","Non-dual interpolated tracks","JEMS log M = " * string(modeltrack1.logM)])
axup = Axis(fig[0,1],xreversed=true, ylabel=L"$\log (L/L_\odot)$") ; axup.xticklabelsvisible = true
zoomin_logM = -0.8
interpolTrack = interpolTrackxs[zoomin_logM]
de.plot!(interpolTrack, axup; color=:lightgreen, scatter=false,linewidth=5)
interpolTrack_dotter = interpolTracks_dotter[zoomin_logM]
de.plot!(interpolTrack_dotter, axup; color=:purple, scatter=false,linewidth=5,linestyle=:dashdot)
text!(axup, interpolTrack.logT_val[1]+0.03,interpolTrack.logL_val[1]+0.05; text = string(zoomin_logM),fontsize=50)
de.plot!(modeltracks[zoomin_logM], axup; color=:black, scatter=false,linewidth=2)
fig
##
figpath = "Figures/DSE_interpolation_vs_dotter.png"; save(figpath, fig, px_per_unit=3); @show figpath
##



################## LINEAR ISOCHRONES
gridpath = "DualRuns/DualGrid3"
get_logM(path) = parse(Float64, split(split(path, "logM_")[2], "_")[1])
function DualGrid(gridpath::String, Xzams_ratio::Float64, Xtams::Float64, nbZeta=1000)
    historypaths = filter(x -> occursin(".history.hdf5", x), readdir(gridpath))
    profilepaths = filter(x -> occursin(".profiles.hdf5", x), readdir(gridpath))
    N = length(historypaths)
    println("$N history files found")
    models = Dict{Float64,de.Model}()
    modeltracks = Dict{Float64, de.Track}()
    X_dual         = ForwardDiff.Dual{}(0.7381,  0.0,1.0,0.0,0.0,0.0)
    Z_dual         = ForwardDiff.Dual{}(0.0134,  0.0,0.0,1.0,0.0,0.0)
    Dfraction_dual = ForwardDiff.Dual{}(0.000312,0.0,0.0,0.0,1.0,0.0)
    R_dual         = ForwardDiff.Dual{}(100*RSUN,0.0,0.0,0.0,0.0,1.0)
    inititial_params_names = [:logM, :X, :Z, :Dfraction, :R]
    all_logMs = []
    only_use_these_logMs = -1.0:0.05:1.0 # [0.0,0.02]; only_use_these_logMs = [0.0]
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
            X_init = Xzams_ratio*model.history_value.X_center[1];
            track = de.Track(model,X_init, Xtams, nbZeta)
        catch 
            println(" Track FAILED for logM = $logM")
            continue
        end
        println(" Track OK for logM = $logM")
        modeltracks[logM] = track
        push!(all_logMs, logM)
    end
    return de.DualGrid(modeltracks, models, all_logMs)
end
##
dualGrid = DualGrid(gridpath, 0.99, 0.01,150);
isochrones = de.Isochrones(dualGrid,15);
isochrones_linear = de.Isochrones(dualGrid,15,do_linear_dotter = true);
##
fig = Figure(size=(1000,800))
ax1 = Axis(fig[1,1], ylabel = L"$\log (L / L_\odot)$",  
xlabel = L"$\log (T / K)$", xreversed = true)

#catter!(ax1, exp10.(isochrones.logT[0.0]), isochrones.logL[0.0])
#catter!(ax1, exp10.(isochrones.logT[1.0]), isochrones.logL[1.0])
#de.plot!(isochrones, ax1 ;scatter=false, nbzeta=20, color=:black,linewidth=10)
de.plot!(isochrones_linear, ax1 ;scatter=false, nbzeta=20, color=:black,linewidth=0.5)

wanted = [0.0,-0.8,0.7,0.8,0.85]; wanted = dualGrid.logMs
for logM in wanted
    track = dualGrid.modelTracks[logM]
    lines!(ax1, track.logT_val, track.logL_val; color=:blue, linewidth=3)
    @show logM
    if logM in -1.0:0.1:1.0 
        text!(ax1, track.logT_val[end], track.logL_val[end];text= string(logM), color=:blue)
    end
end
ax1.xgridvisible = ax1.ygridvisible = true
fig
##
figpath = "Figures/DSE_isochrones_linear.png"; save(figpath, fig, px_per_unit=3); @show figpath