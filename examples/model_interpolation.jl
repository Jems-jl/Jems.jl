using DataFrames
using ForwardDiff
import ForwardDiff.Dual
using Jems.StellarModels
using Interpolations
using HDF5
using Jems.Constants
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
nb = "highres_FULL_5partials"
path = "Jems.jl/DualRuns/"
path = "DualRuns/"
historypath = path * "history_"*string(nb)*".hdf5"
historypath = path * "history_"*string(nb)*".hdf5"
profilespath = path * "profiles_"*string(nb)*".hdf5"
profilespath = path * "profiles_"*string(nb)*".hdf5"
##
#########################################"""
function get_partial_profile_dataframe_from_hdf5(hdf5_filename, value_name, partials_names)
    value_dataframe = StellarModels.get_profile_dataframe_from_hdf5(hdf5_filename, value_name)
    partial_dataframes = [StellarModels.get_profile_dataframe_from_hdf5(hdf5_filename, partial_name) for partial_name in partials_names]
    df_partial = ForwardDiff.Dual.(value_dataframe, partial_dataframes...)
    return df_partial
end
 #doing some bookkeeping stuff
number_of_partials = 5
profile_names = StellarModels.get_profile_names_from_hdf5(profilespath)#all profiles, regular profiles and dual profiles
value_names = [name for name in profile_names if !occursin("partial", name)]
partial_names_unpacked = [name for name in profile_names if occursin("partial", name)]
partial_names = [[partial_name for partial_name in partial_names_unpacked[lo:lo+number_of_partials-1] ] for lo in 1:number_of_partials:(length(partial_names_unpacked))] 


################################################################################## PROFILE OUTPUT
value_names #this list contains the profile names as before, i.e. just the values, nothing special
partial_names #this list contains lists with the corresponding partial names
i = 2
StellarModels.get_profile_dataframe_from_hdf5(profilespath, value_names[i]) #access the ith profile with actual values, as before
bla = get_partial_profile_dataframe_from_hdf5(profilespath, value_names[i], partial_names[i]) #acces the ith profile, but now with Dual numbers, i.e. containg both the values and the partials  
profiles_dual = [get_partial_profile_dataframe_from_hdf5(profilespath, value_name, partial_names) for (value_name,partial_names) in zip(value_names, partial_names)]
profiles_dual # array of all the profiles in hdf5 dual form
profiles_dual[2] #the profile with dual numbers of the second model
#################################################################################

function get_dual_history_dataframe_from_hdf5(hdf5_filename)
    #This function used two functions that were originally defined for profile handling, but they come in handy here
    names = StellarModels.get_profile_names_from_hdf5(hdf5_filename)
    history_value_name = names[1]
    history_partial_names = names[2:end]
    history_value = StellarModels.get_history_dataframe_from_hdf5(hdf5_filename)
    history_partials = [StellarModels.get_profile_dataframe_from_hdf5(hdf5_filename, name) for name in history_partial_names]
    return ForwardDiff.Dual.(history_value, history_partials...)
end
################################################################################# HISTORY OUTPUT
history = StellarModels.get_history_dataframe_from_hdf5(historypath) #as before
history_dual = get_dual_history_dataframe_from_hdf5(historypath)#dataframe with history in dual numbers
(d -> d.partials).(history_dual) #access the values
bla2 = copy(history_dual)
#################################################################################

struct Model
    history
    history_value
    profiles
    profiles_values
    initial_params::Vector{}
    initial_params_names::Vector{}
    initial_params_dict::Dict{}
end

function D_computer(logLs, logTs)
    #distances = sqrt.((logLs.-logLs[1]).^2 .+ (logTs.-logTs[1]).^2)
    distances = zeros(typeof(logLs[1]),length(logLs))
    distances[1] = zero(logLs[1])
    for i in 2:length(logLs)
        delta = sqrt.( (logLs[i]-logLs[i-1] ).^2 + (logTs[i]-logTs[i-1]).^2)
        distances[i] = distances[i-1] + delta
    end
    distances = distances ./distances[end]

    #values = cumsum(distances) .- distances[1]
    @show distances
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
    model.history[!,"logL"] = log10.(model.history[!, "L_surf"])
    model.history[!,"logT"] = log10.(model.history[!, "T_surf"])
    logL_ZAMS = param1_to_param2(ZAMS_X, model.history, "X_center", "logL")
    logT_ZAMS = param1_to_param2(ZAMS_X, model.history, "X_center", "logT")
    logL_TAMS = param1_to_param2(TAMS_X, model.history, "X_center", "logL")
    logT_TAMS = param1_to_param2(TAMS_X, model.history, "X_center", "logT")
    track_history = copy(model.history[ZAMS_index:TAMS_index,:])
    zetas = D_computer(log10.(model.history.L_surf[ZAMS_index:TAMS_index]), log10.(model.history.T_surf[ZAMS_index:TAMS_index]))
    track_history[!,"zeta"] = zetas
    
    track_history[1,"L_surf"] = 10^logL_ZAMS 
    track_history[1,"T_surf"] = 10^logT_ZAMS 
    track_history[end,"L_surf"] = 10^logL_TAMS 
    track_history[end,"T_surf"] = 10^logT_TAMS 
    zetas = LinRange(0,1,nbpoints)
    logL = param1_to_param2.(zetas[2:end-1], Ref(track_history), "zeta", "logL") 
    pushfirst!(logL,logL_ZAMS); push!(logL, logL_TAMS)
    logT = param1_to_param2.(zetas[2:end-1], Ref(track_history), "zeta", "logT")
    pushfirst!(logT,logT_ZAMS); push!(logT, logT_TAMS)
    track_history_value = (dual -> dual.value).(track_history)
    logL_val = (d -> d.value).(logL); logL_partial = (d -> d.partials[1]).(logL)
    logT_val = (d -> d.value).(logT); logT_partial = (d -> d.partials[1]).(logT)
    return Track(model, ZAMS_index, TAMS_index, logL, logL_val, logL_partial, logT, logT_val, logT_partial, zetas, track_history, track_history_value)
end
function plot_track!(track, ax; scatter = true, label = "Track")
    if scatter
        scatter!(ax, track.logT_val, track.logL_val, label=label)
    else 
        lines!(ax, track.logT_val, track.logL_val,  label=label)
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
function ExtrapolTrack(track::Track, delta_logM)
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
function plot!(extrapolGrid::ExtrapolGrid, ax; scatter = true)
    plotfunc = scatter ? scatter! : lines!
    colors_dic = extrapolGrid.colors_dic
    plotfunc(ax, extrapolGrid.track.logT_val, extrapolGrid.track.logL_val, color=:black, label="Original track")
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

function extrapolate(model::Model, delta_params)
    # given a Dual number Dual(value_old, partial_logM, partial_X, partial_Z, ...), we compute in a Taylor way
    # the new value value_new = value_old + [Delta_logM, Delta_X, Delta_Z, ...] * [partial_logM, partial_X, partial_Z, ...]
    # we do this for each value in each profile and for each value in the history
    history_new = copy(model.history)
    history_new =  (dual -> dual.value + sum( delta_params .*dual.partials ) ).(history_new)
    profiles_new = [copy(profile) for profile in model.profiles]
    profiles_new = [(dual -> dual.value + sum( delta_params .*dual.partials ) ).(profile) for profile in profiles_new]
    return Model(nothing, history_new, nothing, profiles_new, model.initial_params, model.initial_params_names, model.initial_params_dict)
end
function star_age_to_X(star_age,history) # computes the X value at a given star age, by linear interpolation
    two_model_numbers = [0,0]
    for modelnr in history.model_number[1].value:history.model_number[end].value
       if modelnr_to_star_age(modelnr) > star_age
            two_model_numbers = [Int(modelnr-1), Int(modelnr)]
            break
       end
    end
    return linear_interpolation(modelnr_to_star_age.(two_model_numbers), modelnr_to_X.(two_model_numbers))(star_age)
end
function X_to_star_age(X,history)
    two_model_numbers = [0,0]
    for modelnr in history.model_number[1].value:history.model_number[end].value
       if modelnr_to_X(modelnr) < X
            two_model_numbers = [Int(modelnr), Int(modelnr-1)]
            break
       end
    end
    return linear_interpolation(modelnr_to_X.(two_model_numbers), modelnr_to_star_age.(two_model_numbers))(X)
end
function find_index(param_value, history, param_name)
    _,index = findmin(abs.(param_value .- history[!,param_name]))
    return index
end
# THE interpolator function
function param1_to_param2(param1_value,history,param1_name,param2_name)
    index_closest = find_index(param1_value, history, param1_name)
    #zetas =  (d->d.value).(history[!,param1_name][1:10])
    if index_closest == 1
        indices = [1,2]
        if history[!,param1_name][2] < history[!, param1_name][1]
            reverse!(indices)
        end
        return linear_interpolation(history[!,param1_name][indices], history[!,param2_name][indices])(param1_value)
    end
    if index_closest == length(history[!,param1_name])
        indices = [length(history[!,param1_name])-1, length(history[!,param1_name])]
        if history[!,param1_name][end] < history[!, param1_name][end-1]
            reverse!(indices)
        end
        return linear_interpolation(history[!,param1_name][indices], history[!,param2_name][indices])(param1_value)
    end
    if history[!,param1_name][index_closest - 1] <= param1_value <= history[!,param1_name][index_closest]
        indices = [index_closest - 1, index_closest]
    elseif history[!, param1_name][index_closest - 1] >= param1_value >= history[!,param1_name][index_closest]
        indices = [index_closest, index_closest - 1]
    elseif history[!,param1_name][index_closest] < param1_value < history[!,param1_name][index_closest + 1]
        indices = [index_closest, index_closest + 1]
    elseif history[!, param1_name][index_closest] > param1_value > history[!,param1_name][index_closest + 1]
        indices = [index_closest + 1, index_closest]
    end
    #@show indices
    #@show history[!,param1_name][indices[1]]
    #@show history[!,param2_name][indices[1]]
    #@show history[!,param1_name][indices[2]]
    #@show history[!,param2_name][indices[2]]
    return my_linear_interpolation(history[!,param1_name][indices], history[!,param2_name][indices])(param1_value)
end

function my_linear_interpolation(xs, ys)
    slope = (ys[2] - ys[1]) / (xs[2] - xs[1])
    return x -> ys[1] + slope * (x - xs[1])
end
##
######### DEFINE MODEL ###################################################################################### DEFINE MODEL
logM_dual      = ForwardDiff.Dual{}(0.0,     1.0,0.0,0.0,0.0,0.0)
mass_dual      = MSUN*10^logM_dual
X_dual         = ForwardDiff.Dual{}(0.7154,  0.0,1.0,0.0,0.0,0.0)
Z_dual         = ForwardDiff.Dual{}(0.0142,  0.0,0.0,1.0,0.0,0.0)
Dfraction_dual = ForwardDiff.Dual{}(0.0,     0.0,0.0,0.0,1.0,0.0)
R_dual         = ForwardDiff.Dual{}(100*RSUN,0.0,0.0,0.0,0.0,1.0)
initial_params = [logM_dual, X_dual, Z_dual, Dfraction_dual, R_dual]
inititial_params_names = [:logM, :X, :Z, :Dfraction, :R]
model1 = Model_constructor(history_dual, profiles_dual, initial_params, inititial_params_names) #construct model!
model1.initial_params_dict[:logM] # access the initial dual number
model1.initial_params_dict[:logM].partials
model1.history # acces history in dual form
model1.history_value # acc es history in value form
model1.profiles # acces list of profiles in dual form
model1.profiles_values # acces list of profiles in value form
delta_params = [1.0,0.0,0.0,0.0,0.0] # delta to use in Taylor expansion
model1_extrapolate = extrapolate(model1, delta_params)
model1_extrapolate.history_value # acces history in value form (extrapolated models DO NOT have history in dual form!)

##
track1 = Track(model1, 0.999*model1.history_value.X_center[1], 0.4*model1.history_value.X_center[1],5000);
#extrapolGrid = ExtrapolGrid(track1, [-0.02,-0.01,0.01,0.02]);
extrapolGrid = ExtrapolGrid(track1, .-[-0.01,-0.02,-0.03,-0.05,-0.08,-0.1,-0.2]);
extrapolGrid = ExtrapolGrid(track1, [-0.01,-0.02,-0.03,-0.05,-0.08,-0.1,-0.2]);
## PLOT TRACKS
fig = Figure();
ax = Axis(fig[1,1], xlabel = L"$\log(Tw_{\text{eff}} / K)$", ylabel = L"$\log (L / L_\odot)$", xreversed = true)
plot!(extrapolGrid,ax)
axislegend(ax,position=:lt)
fig
## PLOT AS FUNCTION OF zeta
fig = Figure(figsize=(2000, 1500))
ax1 = Axis(fig[1,1], ylabel = L"$\log (L / L_\odot)$")
ax2 = Axis(fig[1,2], ylabel = L"$\log (T_{\text{eff}} / K)$")
scatter!(ax1, track1.zeta, track1.logL_val, label = "Original Track", color=:black)
for track in extrapolGrid.extrapoltracks
    scatter!(ax1, track.zeta, track.logL_val, label = "Extrapolated Track", color=extrapolGrid.colors_dic[track])
    scatter!(ax2, track.zeta, track.logT_val, label = "Extrapolated Track", color=extrapolGrid.colors_dic[track])
end

scatter!(ax2, track1.zeta, track1.logT_val, label = "Original Track", color=:black)
ax3 = Axis(fig[2,1], xlabel = L"$\zeta$", ylabel = L"\partial \log L / \partial \log M")
scatter!(ax3, track1.zeta, track1.logL_partial, label = "Original Track", color=:black)
ax4 = Axis(fig[2,2], xlabel = L"$\zeta$", ylabel = L"\partial \log T / \partial \log M")
hlines!(ax4, [0], color=:black, linestyle=:dot,alpha=0.4)
scatter!(ax4, track1.zeta, track1.logT_partial, label = "Original Track", color=:black)
ax5 = Axis(fig[1,3], xlabel=L"$\zeta$", ylabel=L"$X$ and $Y$")
lines!(ax5, track1.history_value.zeta, track1.history_value.X_center, color=:black)
lines!(ax5, track1.history_value.zeta, track1.history_value.Y_center, color=:black,linestyle=:dot)

fig
## PLOT ZAMS
deltas = collect(-0.5:0.01:0.5)
fig = Figure()
ax = Axis(fig[1,1], xlabel = L"$\log(T_{\text{eff}} / K)$", ylabel = L"$\log (L / L_\odot)$", xreversed = true,title="ZAMS and TAMS")
logM = model1.initial_params_dict[:logM].value
plot_track!(track1, ax, label = L"Original track $\log M = %$logM $",scatter=true)

for (i,delta) in enumerate(deltas)
    extrapol = ExtrapolTrack(track1, delta)
    label = L"$\Delta \log M = %$delta $"
    scatter!(ax, extrapol.logT_val[1], extrapol.logL_val[1], label=label,color=:green)
    scatter!(ax, extrapol.logT_val[end], extrapol.logL_val[end], label=label,color=:red)
    #scatter!(ax, extrapol.logT_val[end-1], extrapol.logL_val[end-1], label=label,color=:red)
    if i%10 == 0
        scatter!(ax, extrapol.logT_val, extrapol.logL_val, label=label,color=:gray,alpha=0.4)
    end
end
#axislegend(ax,position=:lb)
fig
##
1+1

###########################################################################################################################
## TESTING the interpolator
# doing some two way checks
param1_to_param2(0.5467236365491246,  model1.history, "X_center", "star_age")
param1_to_param2(5.005870929011321e8,                 model1.history,   "star_age"  , "X_center")
param1_to_param2(0.58,  model1.history, "Y_center", "X_center")
param1_to_param2(0.40557679520000023, model1.history,   "X_center"  , "Y_center")
param1_to_param2(500,  model1.history, "P_surf", "L_surf")
param1_to_param2(25.906216715135315,  model1.history, "L_surf", "P_surf")

# comparing interpolation with closest model
param1_to_param2(0.5,  model1.history, "X_center", "L_surf")  
find_index(0.5, model1.history,"X_center")
                                  model1.history.L_surf[837]  
log10(param1_to_param2(0.7,  model1.history, "X_center", "L_surf"))
find_index(0.7, model1.history,"X_center")
                                  log10(model1.history.L_surf[589])
param1_to_param2(0.5,  model1.history, "X_center", "R_surf")  
                                  model1.history.R_surf[593]  
param1_to_param2(0.5,  model1.history, "X_center", "T_surf")  
                                  model1.history.T_surf[593]  
param1_to_param2(0.5,  model1.history, "Y_center", "X_center")
param1_to_param2(0.5,  model1.history, "X_center", "Y_center")
                                  model1.history.Y_center[593]
param1_to_param2(0.5,  model1.history, "X_center", "X_center")
                                  model1.history.X_center[593]
param1_to_param2(3000,  model1.history, "T_surf", "L_surf")
                                model1.history.L_surf[211] 
param1_to_param2(3000,  model1.history, "T_surf", "T_surf")
                                model1.history.T_surf[211] 
param1_to_param2(18,  model1.history, "L_surf", "L_surf") 
param1_to_param2(500,  model1.history, "P_surf", "P_surf")
##
#plot Lsurf vs X_center0
f = Figure();
ax = Axis(f[1,1])
lines!(ax, model1.history_value.X_center, model1.history_value.L_surf, color = :blue)
f
##
function extrapolate_logM_fixedX(model, delta_logM, X_fixed)
    L_dual_old = param1_to_param2(X_fixed,model.history,"X_center","L_surf")
    logL_dual_old = log10(L_dual_old)
    logL_new = logL_dual_old.value + delta_logM * logL_dual_old.partials[1]
    return 10^logL_new
end

function extrapolate_logM_fixedX_T(model, delta_logM, X_fixed)
    temp_old = param1_to_param2(X_fixed,model.history,"X_center","T_surf")
    temp_new = temp_old.value + delta_logM * temp_old.partials[1]
    return temp_new
end



extrapolate_master(model1, 1, 0.1, "X_center", 0.5, "T_surf")
extrapolate_logM_fixedX_T(model1, 0.1, 0.5)
extrapolate_master_log(model1, 1, 0.1, "X_center", 0.5, "L_surf")
extrapolate_logM_fixedX(model1, 0.1, 0.5)

extrapolate_master_log.(Ref(model1), 1, 0.1, "X_center", [0.5,0.6], "L_surf")
extrapolate_logM_fixedX.(Ref(model1), 0.1, [0.5,0.6])
##
Xarray = collect(0.0001:0.00001:0.9999*model1.history_value.X_center[1])
Xarray = LinRange(0.0001,0.9999*model1.history_value.X_center[1],10000)
Ls = extrapolate_logM_fixedX.(Ref(model1),0.0,Xarray) 
extrapolate_logM_fixedX(model1,0.1,0.9999*model1.history_value.X_center[1])
Temps = extrapolate_logM_fixedX_T.(Ref(model1),0.0,Xarray)
##
f = Figure(); 
ax = Axis(f[1,1])
scatter!(ax, Xarray, Ls, color = :blue)
#scatter!(ax, model1.history_value.T_surf, model1.history_value.L_surf, color = :red)
xlims!(ax,0,0.04)
ylims!(ax,15.1,15.4)
f
##

f = Figure();
ax = Axis(f[1,1])
scatter!(ax, Temps, Ls, color = :blue)
#scatter!(ax, model1.history_value.T_surf, model1.history_value.L_surf, color = :red)
f

##
DeltaLogM = 0.001
f = Figure();
ax = Axis(f[1,1],xreversed=true,xlabel=L"$\log(T_{\text{eff}} / K)$", ylabel=L"$\log (L / L_\odot)$")
begin_index = find_index(0.999*model1.history.X_center[1], model1.history,"X_center")
end_index = find_index(0.01*model1.history.X_center[1], model1.history,"X_center")
scatter!(ax, log10.(model1.history_value.T_surf)[begin_index:end_index], log10.(model1.history_value.L_surf)[begin_index:end_index], color = :blue, label="Original track")
Xarray = collect(0.01*model1.history_value.X_center[1]:0.00001:0.999*model1.history_value.X_center[1])
#Xarray = LinRange(0.0001,0.9999*model1.history_value.X_center[1],10000)
Ls = extrapolate_logM_fixedX.(Ref(model1),DeltaLogM,Xarray) 
Temps = extrapolate_logM_fixedX_T.(Ref(model1),DeltaLogM,Xarray)
scatter!(ax, log10.(Temps), log10.(Ls), color = :red, label=L"Extrapolated track, $\Delta \log M = %$DeltaLogM $")
axislegend(ax, position=:rb)
f
##
deltaLogM_range = [0, 0.1, 0.2]
f = Figure();
ax = Axis(f[1,1],xreversed=true,xlabel=L"$\log(T_{\text{eff}} / K)$", ylabel=L"$\log (L / L_\odot)$")
hi = find_index(0.99*model1.history.X_center[1], model1.history,"X_center")
lo = find_index(0.1*model1.history.X_center[1], model1.history,"X_center")
scatter!(ax, log10.(model1.history_value.T_surf)[hi:lo], log10.(model1.history_value.L_surf)[hi:lo], color = :black, label="Original track")
D_range = LinRange(model1.history.D[hi].value,model1.history.D[lo].value,1000)
function plot_track(ax, DeltaLogM)
    logLs = extrapolate_master_log.(Ref(model1), 1, DeltaLogM, "D", D_range, "L_surf")
    logTs = extrapolate_master_log.(Ref(model1), 1, DeltaLogM, "D", D_range, "T_surf")
    scatter!(ax, logTs, logLs, label=L"$\Delta \log M = %$DeltaLogM $")
end
for delta_logM in deltaLogM_range
    plot_track(ax, delta_logM)
end
axislegend(ax, position=:lc)
f

##

function interpolate(modelA, modelB, logM_wanted)
    hi = find_index(0.999*modelA.history.X_center[1], modelA.history,"X_center")
    lo = find_index(0.01*modelA.history.X_center[1], modelA.history,"X_center")
    D_diff_A = modelA.history.D[hi].value - modelA.history.D[lo].value
    hi = find_index(0.999*modelB.history.X_center[1], modelB.history,"X_center")
    lo = find_index(0.01*modelB.history.X_center[1], modelB.history,"X_center")
    D_diff_B = modelB.history.D[hi].value - modelB.history.D[lo].value
    D_range = LinRange(0, max(D_diff_A,D_diff_B), 1000)
    delta_logM_A = -modelA.initial_params[1]+logM_wanted
    delta_logM_B = -modelB.initial_params[1]+logM_wanted
    logLs_A = extrapolate_master_log.(Ref(modelA), 1, delta_logM_A , "D", D_range, "L_surf")
    logTs_A = extrapolate_master_log.(Ref(modelA), 1, delta_logM_A , "D", D_range, "T_surf")
    logLs_B = extrapolate_master_log.(Ref(modelB), 1, delta_logM_B , "D", D_range, "L_surf")
    logTs_B = extrapolate_master_log.(Ref(modelB), 1, delta_logM_B , "D", D_range, "T_surf")
    logLs = (logLs_A * ( 1 /(abs(delta_logM_A))) + logLs_B * ( 1 /(abs(delta_logM_B))) ) / (1/(abs(delta_logM_A)) + 1/(abs(delta_logM_B)))
    logTs = (logTs_A * ( 1 /(abs(delta_logM_A))) + logTs_B * ( 1 /(abs(delta_logM_B))) ) / (1/(abs(delta_logM_A)) + 1/(abs(delta_logM_B)))
    return logTs, logLs
end  

DeltaLogM = 0.001
f = Figure();
ax = Axis(f[1,1],xreversed=true,xlabel=L"$\log(T_{\text{eff}} / K)$", ylabel=L"$\log (L / L_\odot)$")
hi = find_index(0.999*model1.history.X_center[1], model1.history,"X_center")
lo = find_index(0.001*model1.history.X_center[1], model1.history,"X_center")
scatter!(ax, log10.(model1.history_value.T_surf)[begin_index:end_index], log10.(model1.history_value.L_surf)[begin_index:end_index], color = :blue, label="Original track")
X_range = LinRange(model1.history_value.X_center[lo],model1.history_value.X_center[hi],1000)
logLs = extrapolate_master_log.(Ref(model1), 1, DeltaLogM, "X_center", X_range, "L_surf")
logTs = extrapolate_master_log.(Ref(model1), 1, DeltaLogM, "X_center", X_range, "T_surf")
scatter!(ax, logLs, logTs, color = :red, label=L"Extrapolated track, $\Delta \log M = %$DeltaLogM $")
axislegend(ax, position=:rb)
f

##
function extrapolate_paramB_at_fixed_paramA(model, delta_logM, paramA_fixed, paramA_name,  paramB_name)
    dual_old = param1_to_param2(paramA_fixed,model.history,paramA_name,paramB_name)
    new = dual_old.value + delta_logM * dual_old.partials[1]
end

extrapolate_paramB_at_fixed_paramA(model1,1,0.9999*model1.history_value.X_center[1],"X_center","L_surf")

##
#initial_params = [ForwardDiff.Dual(0.0,1.0)]
#inititial_params_names = [:logM]


##

star_age_to_X(8e8)
X_to_star_age(0.39470445660597037)
##
using Roots

find_zero
##