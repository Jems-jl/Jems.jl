module DualExtrapolation 

using DataFrames
using ForwardDiff
using Jems.StellarModels
using Interpolations#, Dierckx
using HDF5
using DataInterpolations: CubicSpline
using CairoMakie, LaTeXStrings, MathTeXEngine, Makie.Colors, PlotUtils

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

function bookkeeping(historypath, profilespath, number_of_partials)
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
    logM
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
    model.history_value[!,"logL"] = log10.(model.history_value[!, "L_surf"])#adding log 
    model.history[!,"logT"] = log10.(model.history[!, "T_surf"])#adding log 
    model.history_value[!,"logT"] = log10.(model.history_value[!, "T_surf"])#adding log
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
    logM = model.initial_params_dict[:logM].value
    return Track(model, logM,ZAMS_index, TAMS_index, logL, logL_val, logL_partial, logT, logT_val, logT_partial, zetas, track_history, track_history_value)
end
#cubic_interpolator(x,y) = interpolate(x,y,FritschCarlsonMonotonicInterpolation())
#cubic_interpolator(x,y) = Spline1D(x,y)
cubic_interpolator(x,y) = CubicSpline(y,x)

function plot!(track::Track, ax; scatter = true, kwargs...)
    plotfunc = scatter ? scatter! : lines!
    plotfunc(ax, track.logT_val, track.logL_val; kwargs...)
end
struct ExtrapolTrack
    delta_logM
    logM
    original_model::Model
    original_track::Track
    logL_val
    logT_val
    zeta
end
function ExtrapolTrack(track::Track, delta_logM) #EXTRAPOLATION HAPPENS HERE
    logL_new = track.logL_val .+ delta_logM * track.logL_partial
    logT_new = track.logT_val .+ delta_logM * track.logT_partial
    return ExtrapolTrack(delta_logM, track.logM + delta_logM, track.model, track, logL_new, logT_new, track.zeta)
end
function plot!(extrapolTrack, ax; scatter=true, kwargs...)
    plotfunc = scatter ? scatter! : lines!
    plotfunc(ax, extrapolTrack.logT_val, extrapolTrack.logL_val; kwargs...)
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

struct InterpolTrack
    track_down::Track
    extrapoltrack_down::ExtrapolTrack
    track_up::Track
    extrapoltrack_up::ExtrapolTrack
    logM
    logL_val
    logT_val
    zeta
end

function InterpolTrack(track1::Track, track2::Track, logM_wanted)
    logM1 = track1.logM ; logM2 = track2.logM
    if !in_between(logM1, logM_wanted, logM2)
        throw(ErrorException("Wanted logM is not in between the two tracks, interpolation not possible"))
    end
    if length(track1.zeta) != length(track2.zeta)
        throw(ErrorException("Tracks do not have the same zeta sampling, interpolation not possible"))
    end
    track_down = logM1 < logM2 ? track1 : track2; track_up = logM1 < logM2 ? track2 : track1
    delta_logM_down = logM_wanted - track_down.logM; delta_logM_up = logM_wanted - track_up.logM
    weight_up = 1 / abs(delta_logM_down); weight_down = 1 / abs(delta_logM_up)
    total_weight = weight_up + weight_down
    extrapoltrack_from_down = ExtrapolTrack(track_down, delta_logM_down)
    extrapoltrack_from_up = ExtrapolTrack(track_up, delta_logM_up)
    logL_val = (weight_up * extrapoltrack_from_down.logL_val + weight_down * extrapoltrack_from_up.logL_val) / total_weight
    logT_val = (weight_up * extrapoltrack_from_down.logT_val + weight_down * extrapoltrack_from_up.logT_val) / total_weight
    zeta = track_down.zeta
    return InterpolTrack(track_down, extrapoltrack_from_down, track_up, extrapoltrack_from_up, logM_wanted, logL_val, logT_val, zeta)
end

function plot!(interpolTrack::InterpolTrack, ax; scatter=true, kwargs...)
    plotfunc = scatter ? scatter! : lines!
    plotfunc(ax, interpolTrack.logT_val, interpolTrack.logL_val; kwargs...)

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
    _,index = findmin(@. abs(param_value - history[!,param_name]))
    return index
end
function find_index_raw(param_value, array)
    _,index = findmin(@. abs(param_value - array))
    return index
end

in_between(a,x,b) = a <= x <= b || b <= x <= a #helper function

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


end # end of module