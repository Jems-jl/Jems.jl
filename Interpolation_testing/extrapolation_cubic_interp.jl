using DataFrames
using ForwardDiff
import ForwardDiff.Dual
using Jems.StellarModels
using Interpolations
using HDF5
using Jems.Constants
using CairoMakie, LaTeXStrings, MathTeXEngine, Makie.Colors, PlotUtils
#using GLMakie
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
#GLMakie.activate!()

function get_partial_profile_dataframe_from_hdf5(hdf5_filename, value_name, partials_names)
    value_dataframe = StellarModels.get_profile_dataframe_from_hdf5(hdf5_filename, value_name)
    partial_dataframes = [StellarModels.get_profile_dataframe_from_hdf5(hdf5_filename, partial_name) for partial_name in partials_names]
    df_partial = ForwardDiff.Dual.(value_dataframe, partial_dataframes...)
    return df_partial
end
function get_dual_history_dataframe_from_hdf5(hdf5_filename)
    #This function used two functions that were originally defined for profile handling, but they come in handy here
    names = StellarModels.get_profile_names_from_hdf5(hdf5_filename)
    history_value_name = names[1]
    history_partial_names = names[2:end]
    history_value = StellarModels.get_history_dataframe_from_hdf5(hdf5_filename)
    history_partials = [StellarModels.get_profile_dataframe_from_hdf5(hdf5_filename, name) for name in history_partial_names]
    return ForwardDiff.Dual.(history_value, history_partials...)
end

function bookkeeping(historypath, profilespath)
    number_of_partials = 5
    profile_names = StellarModels.get_profile_names_from_hdf5(profilespath)#all profiles, regular profiles and dual profiles
    value_names = [name for name in profile_names if !occursin("partial", name)]
    partial_names_unpacked = [name for name in profile_names if occursin("partial", name)]
    partial_names = [[partial_name for partial_name in partial_names_unpacked[lo:lo+number_of_partials-1] ] for lo in 1:number_of_partials:(length(partial_names_unpacked))] 


    profiles_dual = [get_partial_profile_dataframe_from_hdf5(profilespath, value_name, partial_names) for (value_name,partial_names) in zip(value_names, partial_names)]

    history = StellarModels.get_history_dataframe_from_hdf5(historypath) #as before
    history_dual = get_dual_history_dataframe_from_hdf5(historypath)#dataframe with history in dual numbers

    return history_dual, profiles_dual
end

struct Model
    history
    history_value
    profiles
    profiles_values
    initial_params::Vector{}
    initial_params_names::Vector{}
    initial_params_dict::Dict{}
end

function D_computer_old(logLs, logTs)
    distances = zeros(typeof(logLs[1]),length(logLs))
    distances[1] = zero(logLs[1])
    for i in 2:length(logLs)
        delta = sqrt( (logLs[i]-logLs[i-1] )^2 + (logTs[i]-logTs[i-1])^2 )
        distances[i] = distances[i-1] + delta

    end

    distances = distances ./distances[end]
    return distances
end

function D_computer(logLs, logTs)
    #smaller covered 'distance in L or T space' means a larger weight for each covered distance
    logL_weight = 1/(maximum(logLs) - minimum(logLs)); logT_weight = 1/(maximum(logTs) - minimum(logTs))
    distances = zeros(typeof(logLs[1]),length(logLs))
    distances[1] = zero(logLs[1])
    for i in 2:length(logLs)
        delta = sqrt( (logLs[i]-logLs[i-1] )^2*logL_weight + (logTs[i]-logTs[i-1])^2*logT_weight )
        distances[i] = distances[i-1] + delta
    end
    distances = distances ./distances[end]
    return distances
end

function Model_constructor(history::DataFrame, profiles, initial_params, initial_params_names)
    initial_params_dict = Dict(zip(initial_params_names, initial_params))
    history_value = (dual -> dual.value).(history)
    profiles_values = [(dual -> dual.value).(profile) for profile in profiles]
    Model(history, history_value, profiles, profiles_values, initial_params, initial_params_names,initial_params_dict)
end

struct Track
    model::Model
    ZAMS_index
    TAMS_index
    logL
    logL_val
    logL_partial #partial with respect to logM
    logT
    logT_val
    logT_partial #partial with respect to logM
    zeta #number between 0 and 1, indicating where on the track we are
    history
    history_value
end


function Track(model, ZAMS_X, TAMS_X, nbpoints=1000)
    ZAMS_index = find_index(ZAMS_X, model.history, "X_center")
    TAMS_index = find_index(TAMS_X, model.history, "X_center")
    model.history[!,"logL"] = log10.(model.history[!, "L_surf"])#adding log 
    model.history[!,"logT"] = log10.(model.history[!, "T_surf"])#adding log 
    logL_ZAMS = param1_to_param2(ZAMS_X, model.history, "X_center", "logL") #linear interpolation
    logT_ZAMS = param1_to_param2(ZAMS_X, model.history, "X_center", "logT") #linear interpolation
    logL_TAMS = param1_to_param2(TAMS_X, model.history, "X_center", "logL") #linear interpolation
    logT_TAMS = param1_to_param2(TAMS_X, model.history, "X_center", "logT") #linear interpolation

    track_history = copy(model.history[ZAMS_index:TAMS_index,:])
    track_history[1,"logL"]   = logL_ZAMS 
    track_history[1,"logT"]   = logT_ZAMS 
    track_history[end,"logL"] = logL_TAMS 
    track_history[end,"logT"] = logT_TAMS 
    zetas = D_computer(track_history[!,"logL"], track_history[!,"logT"])
    track_history[!,"zeta"] = zetas
    zetas = collect(LinRange(0,1,nbpoints))
    #zetas = 0:1/(nbpoints-1):1

    interpolator_logL = cubic_interpolator(track_history.zeta, track_history.logL)
    #@show zetas[2:end-1]
    logL = interpolator_logL.(zetas[2:end-1])
    #logL = param1_to_param2.(zetas[2:end-1], Ref(track_history), "zeta", "logL") 
    pushfirst!(logL,logL_ZAMS); push!(logL, logL_TAMS)

    interpolator_logT = cubic_interpolator(track_history.zeta, track_history.logT)
    logT = interpolator_logT.(zetas[2:end-1])
    #logT = param1_to_param2.(zetas[2:end-1], Ref(track_history), "zeta", "logT")
    pushfirst!(logT,logT_ZAMS); push!(logT, logT_TAMS)

    track_history_value = (dual -> dual.value).(track_history)
    logL_val = (d -> d.value).(logL); logL_partial = (d -> d.partials[1]).(logL)
    logT_val = (d -> d.value).(logT); logT_partial = (d -> d.partials[1]).(logT)
    return Track(model, ZAMS_index, TAMS_index, logL, logL_val, logL_partial, logT, logT_val, logT_partial, zetas, track_history, track_history_value)
end
cubic_interpolator(x,y) = interpolate(x,y,FritschCarlsonMonotonicInterpolation())
function plot_track!(track, ax; scatter = true, label = "Track",color=nothing)
    if scatter
        if color == nothing
            scatter!(ax, track.logT_val, track.logL_val, label=label)
        else
            scatter!(ax, track.logT_val, track.logL_val, label=label,color=color)
        end
    else 
        if color == nothing
            lines!(ax, track.logT_val, track.logL_val,  label=label)
        else
            lines!(ax, track.logT_val, track.logL_val,  label=label,color=color)
        end

    end
end
struct ExtrapolTrack
    delta_logM
    original_M
    original_model::Model
    original_track::Track
    logL_val
    logT_val
    zeta
end
function ExtrapolTrack(track::Track, delta_logM) #EXTRAPOLATION HAPPENS HERE
    logL_new = track.logL_val .+ delta_logM * track.logL_partial
    logT_new = track.logT_val .+ delta_logM * track.logT_partial
    return ExtrapolTrack(delta_logM, track.model.initial_params_dict[:logM], track.model, track, logL_new, logT_new, track.zeta)
end

struct ExtrapolGrid
    track::Track
    extrapoltracks::Vector{ExtrapolTrack}
    delta_logM
    colors_dic
end
function ExtrapolGrid(track::Track, deltas)
    extrapoltracks = [ExtrapolTrack(track, delta) for delta in deltas]
    colors_dic = get_colors(extrapoltracks)
    return ExtrapolGrid(track, extrapoltracks, deltas, colors_dic)
end
function get_colors(extrapoltracks)
    colormap = :viridis
    nb_colors = length(extrapoltracks)
    colors = palette(colormap)
    vals = collect(1:nb_colors).*255 ./ nb_colors
    colors_dic = Dict()
    for (i,extrapoltrack) in enumerate(extrapoltracks)
        colors_dic[extrapoltrack] = colors[round(Int,vals[i])]
    end
    return colors_dic
end
function plot!(extrapolGrid::ExtrapolGrid, ax; scatter = true, plot_original = true)
    plotfunc = scatter ? scatter! : lines!
    colors_dic = extrapolGrid.colors_dic
    if plot_original
        plotfunc(ax, extrapolGrid.track.logT_val, extrapolGrid.track.logL_val, color=:black, label="Original track")
    end
    for extrapoltrack in extrapolGrid.extrapoltracks
        delta_logM = extrapoltrack.delta_logM
        label = L"$\Delta \log M = %$delta_logM $"
        plotfunc(ax, extrapoltrack.logT_val, extrapoltrack.logL_val,color=colors_dic[extrapoltrack],label=label)
    end
end

function extrapolate_master(model, init_param_index, init_param_delta, condition_param_name, condition_param_value, target_param_name)
    dual_old = param1_to_param2(condition_param_value,model.history,condition_param_name,target_param_name)
    dual_new = dual_old.value + init_param_delta * dual_old.partials[init_param_index]
    return dual_new
end

function extrapolate_master_log(model, init_param_index, init_param_log_delta, condition_param_name, condition_param_value, target_param_name)
    dual_old = log10(param1_to_param2(condition_param_value,model.history,condition_param_name,target_param_name))
    dual_new = dual_old.value + init_param_log_delta * dual_old.partials[init_param_index]
    return dual_new
end
function find_index(param_value, history, param_name)
    _,index = findmin(abs.(param_value .- history[!,param_name]))
    return index
end
function find_index_raw(param_value, array)
    _,index = findmin(abs.(param_value .- array))
    return index
end



# THE MAIN INTERPOLATOR FUNCTION
function param1_to_param2(param1_value,history,param1_name,param2_name)
    index_closest = find_index(param1_value, history, param1_name)
    #zetas =  (d->d.value).(history[!,param1_name][1:10])
    if index_closest == 1
        indices = [1,2]
        return linear_interpolation(history[!,param1_name][indices], history[!,param2_name][indices])(param1_value)
    end
    if index_closest == length(history[!,param1_name])
        indices = [length(history[!,param1_name])-1, length(history[!,param1_name])]
        return linear_interpolation(history[!,param1_name][indices], history[!,param2_name][indices])(param1_value)
    end
    in_between(a,x,b) = a <= x <= b || b <= x <= a
    if in_between(history[!,param1_name][index_closest - 1],param1_value,history[!,param1_name][index_closest])
        indices = [index_closest - 1, index_closest]
    elseif in_between(history[!,param1_name][index_closest],param1_value,history[!,param1_name][index_closest + 1])
        indices = [index_closest, index_closest + 1]
    else
        throw(ErrorException("Closest value to $param1_name = $param1_value is not crossed in the dataframe."))
    end
    return my_linear_interpolation(history[!,param1_name][indices], history[!,param2_name][indices])(param1_value)
end

function my_linear_interpolation(xs, ys)
    slope = (ys[2] - ys[1]) / (xs[2] - xs[1])
    return x -> ys[1] + slope * (x - xs[1])
end
function arrow_coords(dual_number_x, dual_number_y, deltaM)
    x_begin = dual_number_x.value; y_begin = dual_number_y.value
    return x_begin, dual_number_x.partials[1]*deltaM, y_begin, dual_number_y.partials[1]*deltaM
end
function plot_arrows!(ax,track,zetas,deltaLogM) 
    interpolator_logL = cubic_interpolator(track.zeta, track.logL)
    interpolator_logT = cubic_interpolator(track.zeta, track.logT)
    logLs = interpolator_logL.(zetas); logTs = interpolator_logT.(zetas)
    #logLs = param1_to_param2.(zetas, Ref(track.history), "zeta", "logL")
    #logTs = param1_to_param2.(zetas, Ref(track.history), "zeta", "logT")
    for (logT, logL) in zip(logTs, logLs)
        x, dx, y, dy = arrow_coords(logT, logL, deltaLogM)
        arrows!(ax,[x], [y],  [dx], [dy])
    end
end


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
for i in 1:N
    historypath = joinpath(gridpath, historypaths[i])
    profilepath = joinpath(gridpath, profilepaths[i])
    history_dual, profiles_dual = bookkeeping(historypath, profilepath)
    @assert get_logM(historypath) == get_logM(profilepath)
    logM = get_logM(historypath)
    logM_dual = ForwardDiff.Dual{}(logM, 1.0,0.0,0.0,0.0,0.0)
    initial_params = [logM_dual, X_dual, Z_dual, Dfraction_dual, R_dual]
    model = Model_constructor(history_dual, profiles_dual, initial_params, inititial_params_names)
    models[logM] = model
    track = nothing
    try 
        track = Track(model, 0.05*model.history_value.X_center[1], 0.001, 1000)
    catch 
        println(" Track FAILED at logM = $logM")
        continue
    end
    println(" Track OK at logM = $logM")
    modeltracks[logM] = track
    push!(all_logMs, logM)
end
## SMALL MASS CHANGES
extrapolGrid_00 = ExtrapolGrid(modeltracks[0.0], .- all_logMs);
extrapolGrid_00 = ExtrapolGrid(modeltracks[0.0], [0.01,0.02,-0.01,0.03]);
extrapolGrid_00 = ExtrapolGrid(modeltracks[0.0], [0.01]);
## 
fig = Figure(figsize=(2000,1500))
ax = Axis(fig[1,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$", title=L"Extrapolation from $\log M = 0.0$",xreversed=true)
plot_track!(modeltracks[0.0], ax ; label = L"JEMS $\log M = 0.0$",color=:black,scatter=false)
#plot_track!(modeltracks[0.01], ax; label = L"JEMS $\log M = 0.01$")
#plot_track!(modeltracks[0.02], ax; label = L"JEMS $\log M = 0.02$")
#plot_track!(modeltracks[-0.01], ax; label = L"JEMS $\log M = -0.01$")
#plot_track!(modeltracks[0.03], ax; label = L"JEMS $\log M = 0.03$")
plot!(extrapolGrid_00, ax; plot_original=false,scatter=true)
leg = Legend(fig[1,2],ax)
plot_arrows!(ax,modeltracks[0.0],0:0.05:1,0.01)
fig

## LARGER MASS CHANGES
extrapolGrid_00 = ExtrapolGrid(modeltracks[0.0], [0.1,0.2,0.3]);
fig = Figure(figsize=(2000,1500))
ax = Axis(fig[1,1], xlabel=L"$\log (T_{\text{eff}} / K)$", ylabel=L"$\log (L/L_\odot)$", title=L"Extrapolation from $\log M = 0.0$",xreversed=true)
plot_track!(modeltracks[0.0], ax ; label = L"JEMS $\log M = 0.0$",color=:black,scatter=false)
plot_track!(modeltracks[0.1], ax; label = L"JEMS $\log M = 0.1$", color=:blue)
plot_track!(modeltracks[0.2], ax; label = L"JEMS $\log M = 0.2$", color=:blue)
plot_track!(modeltracks[0.3], ax; label = L"JEMS $\log M = 0.3$", color=:blue)
plot!(extrapolGrid_00, ax; plot_original=false)
leg = Legend(fig[1,2],ax)


fig

## ZETA diagnostic PLOTS
extrapolGrid = ExtrapolGrid(modeltracks[0.0], all_logMs);
logMs_used = [0.01,0.02,0.1]
logMs_used = all_logMs
extrapolGrid = ExtrapolGrid(modeltracks[0.0],logMs_used);
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

