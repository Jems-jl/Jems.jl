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

gridpath = "DualRuns/DualGrid"
gridpath = "DualRuns/DualGrid2"
gridpath = "DualRuns/DualGrid3"
##
path = "DualRuns/DualGrid/logM_-0.1_X_0.7381_.history.hdf5"
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
##
dualGrid = DualGrid(gridpath, 0.99, 0.01,150);
isochrones = de.Isochrones(dualGrid,15);
##
fig = Figure(size=(1000,800))
ax1 = Axis(fig[1,1], ylabel = L"$\log (L / L_\odot)$",  
xlabel = L"$\log (T / K)$", xreversed = true)

#catter!(ax1, exp10.(isochrones.logT[0.0]), isochrones.logL[0.0])
#catter!(ax1, exp10.(isochrones.logT[1.0]), isochrones.logL[1.0])
de.plot!(isochrones, ax1 ;scatter=false, nbzeta=20, color=:black,linewidth=0.5)
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
savepath = "Figures/isochrones_master.png"; save(savepath, fig; px_per_unit=3); @show savepath
##


######################################## CHECK THE GRID DENSITY
struct Difference_with_JEMS
    logM
    diff_logL
    diff_logT
    diff_pyth
    diff_tot
    diff_max
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
    diff_pyth = sqrt.((diff_logL*factor_L).^2 .+ (diff_logT*factor_T).^2); 
    diff_tot = sum(diff_pyth); diff_max = maximum(diff_pyth)
    return Difference_with_JEMS(interpolTrack.logM, diff_logL, diff_logT, diff_pyth, diff_tot, diff_max, interpolTrack.zeta)
end

function max_diff(gridTrack_down::de.Track, gridTrack_middle::de.Track, gridTrack_up::de.Track; nb_tracks_between = 20)
    logM = gridTrack_middle.logM
    interpolTrack = de.InterpolTrack(gridTrack_down, gridTrack_up, logM)
    diff = Difference_with_JEMS(interpolTrack, gridTrack_middle)
    largest_diff = diff.diff_max
    return largest_diff
end

function max_diff(dualGrid::de.DualGrid, leap = 1)
    logMs = sort(dualGrid.logMs)
    max_diffs = Vector{Float64}()
    max_logMs = Vector{Float64}()
    for i in 1+leap:length(logMs)-leap
        logM_down = logMs[i-leap]; logM_middle = logMs[i]; logM_up = logMs[i+leap]
        if round(logM_up - logM_middle,digits=2) != round(logM_middle - logM_down,digits=2)
            println("Skipping logM = $logM_middle")
            continue
        end
        gridTrack_down = dualGrid.modelTracks[logM_down]
        gridTrack_middle = dualGrid.modelTracks[logM_middle]
        gridTrack_up = dualGrid.modelTracks[logM_up]
        largest_diff = max_diff(gridTrack_down, gridTrack_middle, gridTrack_up)
        push!(max_diffs, largest_diff); push!(max_logMs, logM_middle)
        println(" logM = $logM_middle, largest_diff = $largest_diff")
    end
    return max_logMs, max_diffs
end
##


fig = Figure(size=(900,300))
ax = Axis(fig[1,1], xlabel = L"$\log M$", ylabel = L" $\eta$")
#scatter!(ax, max_logMs, log10.(max_diffs))

max_logMs, max_diffs  =  max_diff(dualGrid,6);
scatter!(ax, max_logMs, max_diffs, markersize=20,color=:black, label = L"$\pm 0.3$")

max_logMs, max_diffs  =  max_diff(dualGrid,4);
scatter!(ax, max_logMs, max_diffs, markersize=20,color=:green, label = L"$\pm 0.2$",marker=:cross)

max_logMs, max_diffs  =  max_diff(dualGrid,2);
scatter!(ax, max_logMs, max_diffs, markersize=20,color=:red, label = L"$\pm 0.1$",marker=:x)

max_logMs, max_diffs  =  max_diff(dualGrid,1);
scatter!(ax, max_logMs, max_diffs, markersize=15,color=:blue, label = L"$\pm 0.05$",marker=:diamond)



vlines!(ax, dualGrid.logMs, color=:black, linewidth=2, alpha=0.1)
leg = Legend(fig[1,2], ax)
fig
##
savepath = "Figures/DSE_max_diffs_ifv_mass.png"; save(savepath, fig; px_per_unit=3); @show savepath
##