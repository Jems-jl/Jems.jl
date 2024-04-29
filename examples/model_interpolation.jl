using DataFrames
using ForwardDiff
using Jems.StellarModels
using HDF5
nb = "5partials"
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
    @show names
    history_value_name = names[1]
    history_partial_names = names[2:end]
    history_value = StellarModels.get_history_dataframe_from_hdf5(hdf5_filename)
    history_partials = [StellarModels.get_profile_dataframe_from_hdf5(hdf5_filename, name) for name in history_partial_names]
    return ForwardDiff.Dual.(history_value, history_partials...)
end

################################################################################# HISTORY OUTPUT
history = StellarModels.get_history_dataframe_from_hdf5(historypath) #as before
history_dual = get_dual_history_dataframe_from_hdf5(historypath)#dataframe with history in dual numbers
bla = history_dual
(d -> d.partials).(bla) #access the values
bla2 = copy(bla)
(d ->0.0).(bla2) #access the partials
#################################################################################


struct Model
    history::DataFrame
    profiles::Vector{}
    initial_params::Vector{}
    initial_params_names::Vector{}
    initial_params_dict::Dict{}
end

function Model(history::DataFrame, profiles, initial_params, initial_params_names)
    initial_params_dict = Dict(zip(initial_params_names, initial_params))
    Model(history, profiles, initial_params, initial_params_names, Dict(zip(initial_params_names, initial_params)))
end

function extrapolate(model::Model, delta_params)
    history_new = copy(model.history)
    history_new =  (dual -> dual.value + sum( delta_params .*dual.partials ) ).(history_new)
    profiles_new = [copy(profile) for profile in model.profiles]
    profiles_new = [(dual -> dual.value + sum( delta_params .*dual.partials ) ).(profile) for profile in profiles_new]
    return Model(history_new, profiles_new, model.initial_params, model.initial_params_names)
end
##
#initial_params = [ForwardDiff.Dual(0.0,1.0)]
#inititial_params_names = [:logM]
logM_dual      = ForwardDiff.Dual{tag_external}(0.0,     1.0,0.0,0.0,0.0,0.0)
mass_dual      = MSUN*10^logM_dual
X_dual         = ForwardDiff.Dual{tag_external}(0.7154,  0.0,1.0,0.0,0.0,0.0)
Z_dual         = ForwardDiff.Dual{tag_external}(0.0142,  0.0,0.0,1.0,0.0,0.0)
Dfraction_dual = ForwardDiff.Dual{tag_external}(0.0,     0.0,0.0,0.0,1.0,0.0)
R_dual         = ForwardDiff.Dual{tag_external}(100*RSUN,0.0,0.0,0.0,0.0,1.0)
initial_params = [logM_dual, X_dual, Z_dual, Dfraction_dual, R_dual]
inititial_params_names = [:logM, :X, :Z, :Dfraction, :R]
model1 = Model(history_dual, profiles_dual, initial_params, inititial_params_names)
model1.initial_params_dict[:logM] # access the value
model1.initial_params_dict[:logM].partials
model1.history # acces history in dual form
model1.profiles # acces list of profiles in dual form

delta_params = [1.0,0.0,0.0,0.0,0.0] # delta to use in Taylor expansion
model1_extrapolate = extrapolate(model1, delta_params)
model1_extrapolate.history
(d-> d.value).(bla)
( d-> d.value + sum([5].*d.partials) ).(bla)
