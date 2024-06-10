using LinearAlgebra
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
function cubic_coefficients(y1, partial_y1, y2, partial_y2, Δx)
    A = [Δx^2 Δx^3 ; 2*Δx 3*Δx^2]
    b = [y2 - y1 - Δx*partial_y1 ; partial_y2 - partial_y1]
    x = A \ b; #println("matrix solved")
    return [y1, partial_y1, x[1], x[2]]
end

function cubic_interpolation_function(y1, partial_y1, y2, partial_y2, x1, x2)
    Δx = x2 - x1
    coefficients = cubic_coefficients(y1, partial_y1, y2, partial_y2, Δx)
    function cubic_interpolation(x)
        return coefficients[1] + coefficients[2]*(x-x1) + coefficients[3]*(x-x1)^2 + coefficients[4]*(x-x1)^3
    end
    return cubic_interpolation
end

cubic_coefficients(1, 0, 2, 0, 1)

 
## TEST HET WERKT YES!
fig = Figure();
ax = Axis(fig[1, 1])
random_x = [0.0, 1.0,2.0,3.0]
random_y = [0.0, 18.0,20,10]; @show random_y
random_partial_y = [-0,0,0,0]

cubic_interpolation = cubic_interpolation_function(random_y[1], random_partial_y[1], random_y[2], random_partial_y[2], random_x[1], random_x[2])
x = collect(range(random_x[1], random_x[end],100))
y = [cubic_interpolation(x_i) for x_i in x]
lines!(ax, x, y, color=:blue)
scatter!(ax, random_x, random_y, color=:red, markersize=25)
fig
##

function InterpolTrack_new_cubic(modeltrack1::de.Track, modeltrack2::de.Track, logM_wanted)
    logM1 = modeltrack1.logM ; logM2 = modeltrack2.logM
    if ! de.in_between(logM1, logM_wanted, logM2)
        throw(ErrorException("Wanted logM is not in between the two tracks, interpolation not possible"))
    end
    if length(modeltrack1.zeta) != length(modeltrack2.zeta)
        throw(ErrorException("Tracks do not have the same zeta sampling, interpolation not possible"))
    end
    track_down = logM1 < logM2 ? modeltrack1 : modeltrack2; track_up = logM1 < logM2 ? modeltrack2 : modeltrack1
    delta_logM_down = logM_wanted - track_down.logM; delta_logM_up = logM_wanted - track_up.logM
    logL_vals = Vector{Float64}()
    for (logL_down_val, logL_down_partial, logL_up_val, logL_up_partial) in zip(track_down.logL_val, track_down.logL_partial, track_up.logL_val, track_up.logL_partial)
        cubic_interpolator = cubic_interpolation_function(logL_down_val, logL_down_partial, logL_up_val, logL_up_partial, track_down.logM, track_up.logM)
        logL_val = cubic_interpolator(logM_wanted);
        push!(logL_vals, logL_val)
    end
    logT_vals = Vector{Float64}()
    for (logT_down_val, logT_down_partial, logT_up_val, logT_up_partial) in zip(track_down.logT_val, track_down.logT_partial, track_up.logT_val, track_up.logT_partial)
        cubic_interpolator = cubic_interpolation_function(logT_down_val, logT_down_partial, logT_up_val, logT_up_partial, track_down.logM, track_up.logM)
        logT_val = cubic_interpolator(logM_wanted); push!(logT_vals, logT_val)
    end
    zeta = track_down.zeta
    return de.InterpolTrack(track_down, nothing, track_up, nothing, logM_wanted, logL_vals, logT_vals, zeta)
end


##
interpol_asked = [0.02,0.03,0.04]
interpol_asked = [0.05,0.1,0.15]
interpolTrackxs = Dict()
modeltrack1, modeltrack2 = modeltracks[0.0], modeltracks[0.2];

InterpolTrack_new_cubic(modeltrack1, modeltrack2, 0.1);
for logM in interpol_asked
    interpolTrack = InterpolTrack_new_cubic(modeltrack1, modeltrack2,logM);
    interpolTrackxs[logM] = interpolTrack
end
## EERSTE PLOT INTERPOLATIE IN THESIS
fig = Figure(size=(1500,1200))
ax = Axis(fig[1:2,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$",xreversed=true)
down = de.plot!(modeltrack1, ax; color=:red, scatter=false,linestyle=:dash, linewidth=5)
up = de.plot!(modeltrack2, ax; color=:blue, scatter=false,linestyle=:dot, linewidth=5)
jems_compare = nothing; interpol = nothing; extrapol_down = nothing; extrapol_up = nothing
for logM in interpol_asked
    interpolTrack = interpolTrackxs[logM]
    interpol = de.plot!(interpolTrack, ax; color=:lightgreen, scatter=false,linewidth=5)
    jems_compare = de.plot!(modeltracks[logM], ax; color=:black, scatter=false,linewidth=2)
    text!(ax, interpolTrack.logT_val[1]+0.01,interpolTrack.logL_val[1]-0.01; text = string(logM),fontsize=40)
end

#axislegend(ax,[up, extrapol_up, jems_compare, interpol, extrapol_down, down],["JEMS log M = " * string(modeltrack2.logM),"Extrapolation from above","JEMS tracks for comparison", "Interpolated tracks","Extrapolation from below","JEMS log M = " * string(modeltrack1.logM)])
axup = Axis(fig[0,1],xreversed=true, ylabel=L"$\log (L/L_\odot)$") ; axup.xticklabelsvisible = true
interpolTrack = interpolTrackxs[0.1]
de.plot!(interpolTrack, axup; color=:lightgreen, scatter=false,linewidth=5)
text!(axup, interpolTrack.logT_val[1],interpolTrack.logL_val[1]+0.05; text = string(0.1),fontsize=50)
de.plot!(modeltracks[0.1], axup; color=:black, scatter=false,linewidth=2)
fig
##


## ANDERE MASSA'S #kleurtjes
interpol_asked = 0.05:0.05:0.65; interpolTrackxs = Dict(); diffs_with_JEMS = Dict()
modeltrack1, modeltrack2 = modeltracks[0.0], modeltracks[0.7]
for logM in interpol_asked
    interpolTrack = InterpolTrack_new_cubic(modeltrack1, modeltrack2,logM);
    interpolTrackxs[logM] = interpolTrack; diffs_with_JEMS[logM] = Difference_with_JEMS(interpolTrack, modeltracks[logM])
end
##

fig = Figure(size=(1000,600)); logM1 = modeltrack1.logM; logM2 = modeltrack2.logM; title = "Interpolation between $logM1 and $logM2"
ax = Axis(fig[1,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$",xreversed=true)
de.plot!(modeltrack2, ax; color=:blue, scatter=false,label = "JEMS log M = " * string(modeltrack2.logM),linestyle=:dot, linewidth=5)
scatter!(modeltrack2.logT_val[1], modeltrack2.logL_val[1]; color=:blue, markersize=20)
de.plot!(modeltrack1, ax; color=:red, scatter=false, label = "JEMS log M = " * string(modeltrack1.logM),linestyle=:dash, linewidth=5)
scatter!(modeltrack1.logT_val[1], modeltrack1.logL_val[1]; color=:red, markersize=20)

#lines!(ax, [modeltrack1.logT_val[1],modeltrack2.logT_val[1]], [modeltrack1.logL_val[1],modeltrack2.logL_val[1]]; color=:black, linestyle=:dash,linewidth=2)
zamslogT = Vector{Float64}(); zamslogL = Vector{Float64}(); 
for logM in range(logM1, logM2, 100) #make ZAMS line
    interpolTrack = InterpolTrack_new_cubic(modeltrack1, modeltrack2,logM);
    push!(zamslogT, interpolTrack.logT_val[1]); push!(zamslogL, interpolTrack.logL_val[1])
end
lines!(ax, zamslogT, zamslogL; color=:black, linewidth=2, linestyle=:dash,label="Interpolated ZAMS")


colormap = :viridis; nb_colors = length(interpol_asked); colors = palette(colormap); vals = collect(1:nb_colors).*255 ./ nb_colors
colors = [colors[round(Int,vals[i])] for i in 1:nb_colors]
for (logM,color) in zip(interpol_asked,colors)
    interpolTrack = interpolTrackxs[logM]
    de.plot!(interpolTrack, ax; color=color, scatter=false,linewidth=10,alpha=0.5)
    scatter!(interpolTrack.logT_val[1], interpolTrack.logL_val[1]; color=color, markersize=20)
    de.plot!(modeltracks[logM], ax; color=:black, scatter=false,linewidth=4)
    scatter!(modeltracks[logM].logT_val[1], modeltracks[logM].logL_val[1]; color=:black, markersize=10)
    text!(ax, interpolTrack.logT_val[1]+0.02,interpolTrack.logL_val[1]-0.08; text = string(logM),fontsize=30)
end
ax.xgridvisible = true; ax.ygridvisible = true
#de.plot_arrows!(ax, modeltrack1, 0:0.1:1,0.1);de.plot_arrows!(ax, modeltrack2, 0:0.1:1,-0.1)
axislegend(ax,position = :lb)

#leg = Legend(fig[1,3],ax2;tellwidth=true)
fig
##


## DOTTER Compare
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
    interpolTrack_old = de.InterpolTrack(modeltrack1, modeltrack2,logM);
    interpolTrack = InterpolTrack_new_cubic(modeltrack1, modeltrack2,logM);
    interpolTrack_dotter = de.InterpolTrack_dotter(modeltrack1, modeltrack2,logM);
    interpolTrackxs[logM] = interpolTrack; interpolTracks_dotter[logM] = interpolTrack_dotter
end
## DOTTER COMPARISON THESIS
fig = Figure(size=(1200,1000))
ax = Axis(fig[1:2,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$",xreversed=true)
down = de.plot!(modeltrack1, ax; color=:red, scatter=false,linestyle=:dash, linewidth=5)
up = de.plot!(modeltrack2, ax; color=:blue, scatter=false,linestyle=:dot, linewidth=5)
jems_compare = nothing; interpol = nothing; dotter = nothing; old = nothing
for logM in interpol_asked
    interpolTrack = interpolTrackxs[logM]
    interpol = de.plot!(interpolTrack, ax; color=:green,alpha=0.5, scatter=false,linewidth=10)
    jems_compare = de.plot!(modeltracks[logM], ax; color=:black, scatter=false,linewidth=2)
    text!(ax, interpolTrack.logT_val[1]+0.04,interpolTrack.logL_val[1]-0.03; text = string(logM),fontsize=40)
    interpolTrack_dotter = interpolTracks_dotter[logM]
    dotter = de.plot!(interpolTrack_dotter, ax; color=:purple, scatter=false,linewidth=5,linestyle=:dashdot)
    old = de.plot!(de.InterpolTrack(modeltrack1, modeltrack2,logM), ax; color=:orange, scatter=false,linewidth=5,linestyle=:dot)
end

axislegend(ax,[up, jems_compare, interpol, old, dotter, down],
                    ["JEMS log M = " * string(modeltrack2.logM),"JEMS tracks for comparison", 
                     "Interpolated tracks (new method)","Interpolated tracks (old method)", 
                     "Non-dual interpolated tracks","JEMS log M = " * string(modeltrack1.logM)])
axup = Axis(fig[0,1],xreversed=true, ylabel=L"$\log (L/L_\odot)$") ; axup.xticklabelsvisible = true
zoomin_logM = -0.8
interpolTrack = interpolTrackxs[zoomin_logM]
de.plot!(interpolTrack, axup; color=:green, scatter=false,linewidth=10,alpha=0.5)
interpolTrack_dotter = interpolTracks_dotter[zoomin_logM]
de.plot!(interpolTrack_dotter, axup; color=:purple, scatter=false,linewidth=5,linestyle=:dashdot)
text!(axup, interpolTrack.logT_val[1]+0.03,interpolTrack.logL_val[1]+0.05; text = string(zoomin_logM),fontsize=50)
de.plot!(modeltracks[zoomin_logM], axup; color=:black, scatter=false,linewidth=2)
de.plot!(de.InterpolTrack(modeltrack1, modeltrack2,zoomin_logM), axup; color=:orange, scatter=false,linewidth=5,linestyle=:dot)
fig

## ########################################
function make_interpoltracks(track1::de.Track, track2::de.Track, nb_masses_between::Int; do_linear_dotter = false)
    logM1 = track1.logM; logM2 = track2.logM
    track_down = logM1 < logM2 ? track1 : track2; track_up = logM1 < logM2 ? track2 : track1
    masses = range(logM1, logM2, length=nb_masses_between+2)[2:end-1] #tja de randen doen nu ook mee
    interpolTracks = Vector{de.InterpolTrack}()
    constructor = InterpolTrack_new_cubic
    for logM in masses
        push!(interpolTracks, constructor(track_down, track_up, logM))
    end
    return interpolTracks
end

function make_interpoltracks(dualGrid::de.DualGrid, nb_masses_between; do_linear_dotter = false)
    logMs = sort(dualGrid.logMs)
    all_interpoltracks = []
    for i in 1:length(logMs)-1
        logM_down = logMs[i]; logM_up = logMs[i+1]
        track1 = dualGrid.modelTracks[logM_down]; track2 = dualGrid.modelTracks[logM_up]
        interpolTracks = make_interpoltracks(track1, track2, nb_masses_between; do_linear_dotter = do_linear_dotter)
        all_interpoltracks = vcat(all_interpoltracks, interpolTracks)
    end
    return all_interpoltracks
end

function Isochrones_new(dualGrid::de.DualGrid, nb_masses_between::Int; massrange = nothing, do_linear_dotter = false)
    if massrange != nothing
        dualGrid = DualGrid_cut(dualGrid, massrange[1], massrange[2])
    end
    all_interpoltracks = make_interpoltracks(dualGrid, nb_masses_between; do_linear_dotter = do_linear_dotter)
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
    return de.Isochrones(zetas, logLs, logTs)
end

##
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
dualGrid = DualGrid(gridpath, 0.99, 0.01,150);
isochrones = de.Isochrones(dualGrid,15);
isochrones_linear = de.Isochrones(dualGrid,15,do_linear_dotter = true);
isochrones_new = Isochrones_new(dualGrid, 15)
##
fig = Figure(size=(1000,800))
ax1 = Axis(fig[1,1], ylabel = L"$\log (L / L_\odot)$",  
xlabel = L"$\log (T / K)$", xreversed = true)

#catter!(ax1, exp10.(isochrones.logT[0.0]), isochrones.logL[0.0])
#catter!(ax1, exp10.(isochrones.logT[1.0]), isochrones.logL[1.0])
#de.plot!(isochrones, ax1 ;scatter=false, nbzeta=20, color=:black,linewidth=10)
#de.plot!(isochrones_linear, ax1 ;scatter=false, nbzeta=20, color=:black,linewidth=0.5)
de.plot!(isochrones_new, ax1 ;scatter=false, nbzeta=20, color=:black,linewidth=0.5)
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
