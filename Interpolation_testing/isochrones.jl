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


struct DualGrid
    modelTracks::Dict{Float64, de.Track}
    models::Dict{Float64, de.Model}
    logMs::Vector{Float64}
end

function DualGrid(gridpath::String, Xzams_ratio::Float64, Xtams::Float64, nbZeta=1000)
    historypaths = filter(x -> occursin(".history.hdf5", x), readdir(gridpath))
    profilepaths = filter(x -> occursin(".profiles.hdf5", x), readdir(gridpath))
    N = length(historypaths)
    println("$N history files found")
    models = Dict() #Dict{Float64,de.Model}()
    modeltracks = Dict() #Dict{Float64, de.Track}()
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
    @show typeof(models)
    @show typeof(modeltracks)
    #return modeltracks, models, all_logMs
    return DualGrid(modeltracks, models, all_logMs)
end


function DualGrid_cut(dualGrid::DualGrid, logM_min::Float64, logM_max::Float64)
    modelTracks = Dict{Float64, de.Track}()
    models = Dict{Float64, de.Model}()
    logMs = []
    for logM in dualGrid.logMs
        if logM_min <= logM <= logM_max
            modelTracks[logM] = dualGrid.modelTracks[logM]
            models[logM] = dualGrid.models[logM]
            push!(logMs, logM)
        end
    end
    return DualGrid(modelTracks, models, logMs)
end

##
dualGrid = DualGrid(gridpath, 0.99, 0.01,5);
##

struct Isochrones
    zeta::Vector{Float64}
    logL::Dict{Float64, Vector{Float64}}
    logT::Dict{Float64, Vector{Float64}}
end

function make_interpoltracks(track1::de.Track, track2::de.Track, nb_masses_between::Int)
    logM1 = track1.logM; logM2 = track2.logM
    track_down = logM1 < logM2 ? track1 : track2; track_up = logM1 < logM2 ? track2 : track1
    masses = range(logM1, logM2, length=nb_masses_between+2)[2:end-1] #tja de randen doen nu ook mee
    interpolTracks = Vector{de.InterpolTrack}()
    for logM in masses
        push!(interpolTracks, de.InterpolTrack(track_down, track_up, logM))
    end
    return interpolTracks
end

function make_interpoltracks(dualGrid::DualGrid, nb_masses_between)
    logMs = sort(dualGrid.logMs)
    all_interpoltracks = []
    for i in 1:length(logMs)-1
        logM_down = logMs[i]; logM_up = logMs[i+1]
        track1 = dualGrid.modelTracks[logM_down]; track2 = dualGrid.modelTracks[logM_up]
        interpolTracks = make_interpoltracks(track1, track2, nb_masses_between)
        all_interpoltracks = vcat(all_interpoltracks, interpolTracks)
    end
    return all_interpoltracks
end

function Isochrones(dualGrid::DualGrid, nb_masses_between::Int; massrange = nothing)
    if massrange != nothing
        dualGrid = DualGrid_cut(dualGrid, massrange[1], massrange[2])
    end
    all_interpoltracks = make_interpoltracks(dualGrid, nb_masses_between)
    logLs = Dict{Float64, Vector{Float64}}(); logTs = Dict{Float64, Vector{Float64}}()
    zetas = all_interpoltracks[1].zeta

    for zeta in zetas #initialize empty arrays
        logLs[zeta] = Vector{Float64}(); logTs[zeta] = Vector{Float64}()
    end
    #loop over all interpolated tracks
    for interpoltrack in all_interpoltracks
        #for each interpolated track, push all logL and logT values in the corresponding zeta
        for (zeta, logL, logT) in zip(interpoltrack.zeta, interpoltrack.logL_val, interpoltrack.logT_val)
            push!(logLs[zeta], logL); push!(logTs[zeta], logT)
        end
    end
    return Isochrones(zetas, logLs, logTs)
end

function plott!(ax, isochrones::Isochrones; scatter = true, kwargs...)
    plotfunc = scatter ? scatter! : lines!
    for zeta in isochrones.zeta
        plotfunc(ax, isochrones.logT[zeta], isochrones.logL[zeta]; kwargs...)
        #plotfunc(ax, exp10.(isochrones.logT[zeta]), exp10.(isochrones.logL[zeta]); kwargs...)
    end
end
##
isochrones = Isochrones(dualGrid,50)
isochrones = Isochrones(dualGrid,5; massrange = [-0.5,0.0])
##
fig = Figure()
ax1 = Axis(fig[1,1], ylabel = L"$\log (L / L_\odot)$",  
xlabel = L"$\log (T / K)$", xreversed = true)

#catter!(ax1, exp10.(isochrones.logT[0.0]), isochrones.logL[0.0])
#catter!(ax1, exp10.(isochrones.logT[1.0]), isochrones.logL[1.0])
plott!(ax1, isochrones;scatter=false,color=:black,linewidth=0.1)
wanted = [0.0,-0.8,0.7,0.8,0.85]; wanted = dualGrid.logMs
for logM in wanted
    track = dualGrid.modelTracks[logM]
    lines!(ax1, track.logT_val, track.logL_val; color=:black, linewidth=5)
end
fig
##