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
    only_use_these_logMs =[]# [0.0,0.02]; only_use_these_logMs = [0.0]
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
isochrones = de.Isochrones(dualGrid,5; massrange = [0.0,0.2])
isochrones = de.Isochrones(dualGrid,15)
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