using DataFrames
using ForwardDiff
using Jems.StellarModels
using Interpolations
using HDF5
using Jems.Constants
using CairoMakie, LaTeXStrings, MathTeXEngine
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
nb = "FULL_5partials"
path = "DualRuns/"
historypath = path * "history_"*string(nb)*".hdf5"
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

function Model_constructor(history::DataFrame, profiles, initial_params, initial_params_names)
    initial_params_dict = Dict(zip(initial_params_names, initial_params))
    history_value = (dual -> dual.value).(history)
    profiles_values = [(dual -> dual.value).(profile) for profile in profiles]
    Model(history, history_value, profiles, profiles_values, initial_params, initial_params_names,initial_params_dict )
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

function extrapolate_logM_fixedX(model, delta_logM, X_fixed)
    param1_to_param2(param1_value,history,param1_name,param2_name)
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
            @show two_model_numbers
            break
       end
    end
    return linear_interpolation(modelnr_to_X.(two_model_numbers), modelnr_to_star_age.(two_model_numbers))(X)
end

# the interpolator function
function param1_to_param2(param1_value,history,param1_name,param2_name)
    _,index = findmin(abs.(param1_value .- history[!,param1_name]))
    indices = [index-1,index+1]
    if history[!,param1_name][indices[2]] < history[!,param1_name][indices[1]]
        reverse!(indices) # reverse order, to keep the lienar_interpolation function happy
    end
    return linear_interpolation(history[!,param1_name][indices], history[!,param2_name][indices])(param1_value)
end
##

param1_to_param2(0.5467236365491246,  model1.history, "X_center", "star_age")
param1_to_param2(5e8,                 model1.history,   "star_age"  , "X_center")
param1_to_param2(0.58,  model1.history, "Y_center", "X_center")
param1_to_param2(0.40580000000000016, model1.history,   "X_center"  , "Y_center")

param1_to_param2(0.5,  model1.history, "X_center", "L_surf")  
                                  model1.history.L_surf[593]
log10(param1_to_param2(0.5,  model1.history, "X_center", "L_surf"))
                                  log10(model1.history.L_surf[593])
param1_to_param2(0.5,  model1.history, "X_center", "R_surf")  
                                  model1.history.R_surf[593]  
param1_to_param2(0.5,  model1.history, "X_center", "T_surf")  
                                  model1.history.T_surf[593]  
param1_to_param2(0.5,  model1.history, "X_center", "Y_center")
                                  model1.history.Y_center[593]
param1_to_param2(3000,  model1.history, "T_surf", "L_surf")   
                                 model1.history.L_surf[211]   
##
#plot Lsurf vs X_center
f = Figure();
ax = Axis(f[1,1])
lines!(ax, model1.history_value.X_center, model1.history_value.L_surf, color = :blue)
f
##
#initial_params = [ForwardDiff.Dual(0.0,1.0)]
#inititial_params_names = [:logM]
###################### DEFINE MODEL
logM_dual      = ForwardDiff.Dual{}(0.0,     1.0,0.0,0.0,0.0,0.0)
mass_dual      = MSUN*10^logM_dual
X_dual         = ForwardDiff.Dual{}(0.7154,  0.0,1.0,0.0,0.0,0.0)
Z_dual         = ForwardDiff.Dual{}(0.0142,  0.0,0.0,1.0,0.0,0.0)
Dfraction_dual = ForwardDiff.Dual{}(0.0,     0.0,0.0,0.0,1.0,0.0)
R_dual         = ForwardDiff.Dual{}(100*RSUN,0.0,0.0,0.0,0.0,1.0)
initial_params = [logM_dual, X_dual, Z_dual, Dfraction_dual, R_dual]
inititial_params_names = [:logM, :X, :Z, :Dfraction, :R]
model1 = Model_constructor(history_dual, profiles_dual, initial_params, inititial_params_names)
model1.initial_params_dict[:logM] # access the initial dual number
model1.initial_params_dict[:logM].partials
model1.history # acces history in dual form
model1.history_value # acces history in value form
model1.profiles # acces list of profiles in dual form
model1.profiles_values # acces list of profiles in value form
delta_params = [1.0,0.0,0.0,0.0,0.0] # delta to use in Taylor expansion
model1_extrapolate = extrapolate(model1, delta_params)
model1_extrapolate.history_value # acces history in value form (extrapolated models DO NOT have history in dual form!)



history = model1.history
modelnr_to_X = modelnr -> history.X_center[Int(modelnr)]
modelnr_to_X(5)
modelnr_to_star_age = modelnr -> history.star_age[Int(modelnr)]
#star_age_to_modelnr = star_age -> floor(Int, ( find_zero(modelnr -> modelnr_to_star_age(modelnr)-star_age,(1,history.model_number[end].value), Bisection() ) ))
#star_age_to_modelnr(1647.0)

##

star_age_to_X(8e8)
X_to_star_age(0.39470445660597037)
##
using Roots

find_zero
##
blaa = modelnr -> modelnr_to_star_age(modelnr)-history.star_age[5]
blaa(4)
##

