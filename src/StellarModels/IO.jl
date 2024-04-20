using HDF5
using DataFrames
using Printf

export add_history_option, add_profile_option,
    history_output_units, history_output_functions,
    profile_output_units, profile_output_functions

const width = 9
const decimals = 4
const floatstr = "%#$width.$decimals" * "g "
const intstr = "%$width" * "i "
const line1fmt = Printf.Format(intstr * floatstr^6 * intstr * "\n")
const line2fmt = Printf.Format(floatstr^7 * intstr * "\n")
const header = """
    model     logdt      logL   logTeff     logPs     logρs    H_cntr     iters
     mass       age      logR     logTc     logPc     logρc   He_cntr     zones
-------------------------------------------------------------------------------
"""

const history_output_units::Dict{String,String} = Dict()
const history_output_functions::Dict{String,Function} = Dict()
const profile_output_units::Dict{String,String} = Dict()
const profile_output_functions::Dict{String,Function} = Dict()

function add_history_option(name, unit, func)
    if haskey(history_output_units, name)
        throw(ArgumentError("Key $name is already part of the history output options"))
    end
    history_output_units[name] = unit
    history_output_functions[name] = func
end

#function get_history_output_function_from_dual(func,sm)
#    prinln("STARTED get_history_output_function_from_dual")
#    dual_variable = func(sm)
#    @show dual_variable
#    @show typeof(dual_variable)
#    return dual_variable

# general properties
add_history_option("star_age", "year", sm -> sm.props.time / SECYEAR)
add_history_option("dt", "year", sm -> sm.props.dt / SECYEAR)
add_history_option("model_number", "unitless", sm -> sm.props.model_number)
add_history_option("star_mass", "Msun", sm -> sm.props.mstar / MSUN)

# surface properties
add_history_option("R_surf", "Rsun", sm -> exp(get_value(sm.props.lnr[sm.props.nz])) / RSUN)
add_history_option("L_surf", "Lsun", sm -> get_value(sm.props.L[sm.props.nz]))
add_history_option("T_surf", "K", sm -> exp(get_value(sm.props.lnT[sm.props.nz])))
add_history_option("ρ_surf", "g*cm^-3", sm -> exp(get_value(sm.props.lnρ[sm.props.nz])))
add_history_option("P_surf", "dyne", sm -> get_value(sm.props.eos_res[sm.props.nz].P))
add_history_option("X_surf", "unitless", sm -> get_value(sm.props.xa[sm.props.nz, sm.network.xa_index[:H1]]))
add_history_option("Y_surf", "unitless", sm -> get_value(sm.props.xa[sm.props.nz, sm.network.xa_index[:He4]]))

# central properties
add_history_option("T_center", "K", sm -> exp(get_value(sm.props.lnT[1])))
add_history_option("ρ_center", "g*cm^-3", sm -> exp(get_value(sm.props.lnρ[1])))
add_history_option("P_center", "dyne", sm -> exp(get_value(sm.props.eos_res[1].P)))
add_history_option("X_center", "unitless", sm -> get_value(sm.props.xa[1, sm.network.xa_index[:H1]]))
add_history_option("Y_center", "unitless", sm -> get_value(sm.props.xa[1, sm.network.xa_index[:He4]]))

function add_profile_option(name, unit, func)
    if haskey(profile_output_units, name)
        throw(ArgumentError("Key $name is already part of the history output options"))
    end
    profile_output_units[name] = unit
    profile_output_functions[name] = func
end

# mesh properties
add_profile_option("zone", "unitless", (sm, k) -> k)
add_profile_option("mass", "Msun", (sm, k) -> sm.props.m[k] / MSUN)
add_profile_option("dm", "Msun", (sm, k) -> sm.props.dm[k] / MSUN)

# thermodynamic properties
add_profile_option("log10_r", "log10(Rsun)", (sm, k) -> get_value(sm.props.lnr[k]) * log10_e - log10(RSUN))
add_profile_option("log10_P", "log10(dyne)", (sm, k) -> get_value(sm.props.eos_res[k].lnP) * log10_e)
add_profile_option("log10_T", "log10(K)", (sm, k) -> get_value(sm.props.lnT[k]) * log10_e)
add_profile_option("log10_ρ", "log10_(g*cm^-3)", (sm, k) -> get_value(sm.props.lnρ[k]) * log10_e)
add_profile_option("luminosity", "Lsun", (sm, k) -> get_value(sm.props.L[k]) / LSUN)

# abundances
add_profile_option("X", "unitless", (sm, k) -> get_value(sm.props.xa[k, sm.network.xa_index[:H1]]))
add_profile_option("Y", "unitless", (sm, k) -> get_value(sm.props.xa[k, sm.network.xa_index[:He4]]))

# temperature gradients
add_profile_option("nabla_a_face", "unitless", (sm, k) -> get_value(sm.props.∇ₐ_face[k]))
add_profile_option("nabla_r_face", "unitless", (sm, k) -> get_value(sm.props.∇ᵣ_face[k]))
add_profile_option("nabla_face", "unitless", (sm, k) -> get_value(sm.props.turb_res[k].∇))
add_profile_option("D_face", "unitless", (sm, k) -> get_value(sm.props.turb_res[k].D_turb))

"""
    create_output_files(sm::StellarModel)

Creates output files for history and profile data
"""
function create_output_files!(sm::StellarModel{TNUMBER, TDUALFULL, TPROPS,
    TEOS,TKAP,TNET, TTURB, TSOLVER}) where {TNUMBER, TDUALFULL, TPROPS, TEOS, TKAP, TNET, TTURB, TSOLVER}
    # Create history file
    sm.history_file = h5open(sm.opt.io.hdf5_history_filename, "w")
    data_cols = sm.opt.io.history_values
    ncols = length(data_cols)

    # verify validity of column names
    for i in eachindex(data_cols)
        if data_cols[i] ∉ keys(history_output_functions)
            throw(ArgumentError("Invalid name for history data column, 
                :$(data_cols[i])"))
        end
    end

    # Create history dataset in HDF5 file
    # Dataset is created with size (0, ncols), we will add rows by using the HDF5.set_extent_dims function
    # the (-1, ncols) is used to define the maximum extent of the dataset, -1 indicates that it is unbound
    # in number of rows. The chunk size is used for compression. Smaller chunk sizes will result in worse
    # compression but faster writes.
    # The compression level can be anywhere between 0 and 9, 0 being no compression 9 being the highest.
    # Compression is lossless.
    #history = create_dataset(sm.history_file, "history", Float64, ((0, ncols), (-1, ncols)),
    #                            chunk=(sm.opt.io.hdf5_history_chunk_size, ncols),
    #                            compress=sm.opt.io.hdf5_history_compression_level)

    #we want a history hdf5 file which will essentially contain several of the above defined tables 
    history = create_dataset(sm.history_file, "history", Float64, ((0, ncols), (-1, ncols)),
                                chunk=(sm.opt.io.hdf5_history_chunk_size, ncols),
                                compress=sm.opt.io.hdf5_history_compression_level)
    
    # next up, include the units for all quantities. No need to recheck columns.
    attrs(history)["column_units"] = [history_output_units[data_cols[i]] for i in eachindex(data_cols)]
    # Finally, place column names
    attrs(history)["column_names"] = [data_cols[i] for i in eachindex(data_cols)]

    if TNUMBER != Float64
        println("TNUMBER is not Float64")
        number_of_partials = TNUMBER.parameters[3]
        dual_histories = [create_dataset(sm.history_file, "dualhistory_$i", Float64, ((0, ncols), (-1, ncols)),
                                chunk=(sm.opt.io.hdf5_history_chunk_size, ncols),
                                compress=sm.opt.io.hdf5_history_compression_level) for i in 1:number_of_partials]
        for dual_history in dual_histories
            attrs(dual_history)["column_units"] = [history_output_units[data_cols[i]] for i in eachindex(data_cols)]
            attrs(dual_history)["column_names"] = [data_cols[i] for i in eachindex(data_cols)]
        end
    end


    if (!sm.opt.io.hdf5_history_keep_open)
        close(sm.history_file)
    end
    

    # Create profile file
    sm.profiles_file = h5open(sm.opt.io.hdf5_profile_filename, "w")
    data_cols = sm.opt.io.profile_values
    # verify validity of column names
    for i in eachindex(data_cols)
        if data_cols[i] ∉ keys(profile_output_functions)
            throw(ArgumentError("Invalid name for history data column,
                :$(data_cols[i])"))
        end
    end
    if (!sm.opt.io.hdf5_profile_keep_open)
        close(sm.profiles_file)
    end
end

function close_output_files!(sm)
    if(sm.opt.io.hdf5_history_keep_open)
        close(sm.history_file)
    end
    if(sm.opt.io.hdf5_profile_keep_open)
        close(sm.profiles_file)
    end
end

"""
    write_data(sm::StellarModel)

Saves data (history/profile) for the current model, as required by the settings in `sm.opt.io`.
"""
function write_data(sm::StellarModel{TNUMBER, TDUALFULL, TPROPS,
    TEOS,TKAP,TNET, TTURB, TSOLVER}) where {TNUMBER, TDUALFULL, TPROPS, TEOS, TKAP, TNET, TTURB, TSOLVER}
    type_in_symbol_form =  TNUMBER.name.name
    # do history
    if (sm.opt.io.history_interval > 0)
        file_exists = isfile(sm.opt.io.hdf5_history_filename)
        if !file_exists  # create file if it doesn't exist yet
            throw(ErrorException("History file does not exist at $(sm.opt.io.hdf5_history_filename)"))
        end
        if (sm.props.model_number % sm.opt.io.history_interval == 0)
            if (!sm.opt.io.hdf5_history_keep_open)
                sm.history_file = h5open(sm.opt.io.hdf5_history_filename, "r+")
            end
            data_cols = sm.opt.io.history_values
            ncols = length(data_cols)

            # after being sure the header is there, print the data
            history = sm.history_file["history"]
            HDF5.set_extent_dims(history, (size(history)[1] + 1, ncols))
            if TNUMBER != Float64
                dual_histories = [sm.history_file["dualhistory_$i"] for i in 1:TNUMBER.parameters[3] ]
                for dual_history in dual_histories
                    HDF5.set_extent_dims(dual_history, (size(dual_history)[1] + 1, ncols))
                end
            end
            for i in eachindex(data_cols)
                colname = data_cols[i]
                if TNUMBER != Float64
                    if colname == "model_number"
                        history[end, i] = history_output_functions[data_cols[i]](sm) #model number is never a dual number
                        for (k,dual_history) in enumerate(dual_histories)
                            dual_history[end,i] = history_output_functions[data_cols[i]](sm)#add model number
                        end
                    else 
                        history[end, i] = history_output_functions[data_cols[i]](sm).value #value for dual numbers
                        for (k,dual_history) in enumerate(dual_histories)
                            dual_history[end,i] = history_output_functions[data_cols[i]](sm).partials[k] #partials for dual numbers
                        end
                    end
                else
                    history[end, i] = history_output_functions[data_cols[i]](sm) #for non-dual (normal) numbers
                end
            end
            if (!sm.opt.io.hdf5_history_keep_open)
                close(sm.history_file)
            end
        end
    end
    # history done
    # do profile
    if (sm.opt.io.profile_interval > 0)
        file_exists = isfile(sm.opt.io.hdf5_profile_filename)
        if !file_exists  # create file if it doesn't exist yet
            throw(ErrorException("Profile file does not exist at $(sm.opt.io.hdf5_profile_filename)"))
        end
        if (sm.props.model_number % sm.opt.io.profile_interval == 0)
            if (!sm.opt.io.hdf5_profile_keep_open)
                sm.profiles_file = h5open(sm.opt.io.hdf5_profile_filename, "r+")
            end
            data_cols = sm.opt.io.profile_values
            ncols = length(data_cols)
            # Save current profile
            profile = create_dataset(sm.profiles_file,
                                        "$(lpad(sm.props.model_number,sm.opt.io.hdf5_profile_dataset_name_zero_padding,"0"))",
                                        Float64, ((sm.props.nz, ncols), (sm.props.nz, ncols));
                                        chunk=(sm.opt.io.hdf5_profile_chunk_size, ncols),
                                        compress=sm.opt.io.hdf5_profile_compression_level)

            # next up, include the units for all quantities. No need to recheck columns.
            attrs(profile)["column_units"] = [profile_output_units[data_cols[i]] for i in eachindex(data_cols)]
            # Place column names
            attrs(profile)["column_names"] = [data_cols[i] for i in eachindex(data_cols)]

            if TNUMBER == Float64
                # store data
                for i in eachindex(data_cols), k = 1:(sm.props.nz)
                    profile[k, i] = profile_output_functions[data_cols[i]](sm, k)
                end
            else TNUMBER != Float64 #then it is a dual number
                for i in eachindex(data_cols), k = 1:(sm.props.nz)
                    #storing the actual profile data, not yet the partials
                    colname = data_cols[i]
                    #@show data_cols[i]
                    #println(" ")
                    #@show profile_output_functions[data_cols[i]](sm, k)
                    #println(" ")
                    #@show typeof(profile_output_functions[data_cols[i]](sm, k))
                    if colname == "zone" #zone is not a dual number
                        profile[k, i] = profile_output_functions[colname](sm, k)
                    else
                        profile[k, i] = profile_output_functions[colname](sm, k).value
                    end
                end
                number_of_partials = TNUMBER.parameters[3]
                dual_profiles = []
                for partial_index in 1:number_of_partials #loop over all input parameters
                    #create a new profile dual_profile for the (partial_index)th partial
                    dual_profile = create_dataset(sm.profiles_file,
                        "$(lpad(sm.props.model_number,sm.opt.io.hdf5_profile_dataset_name_zero_padding,"0"))dual_$partial_index",
                        Float64, ((sm.props.nz, ncols), (sm.props.nz, ncols));
                        chunk=(sm.opt.io.hdf5_profile_chunk_size, ncols),
                        compress=sm.opt.io.hdf5_profile_compression_level)
                    # next up, include the units for all quantities. No need to recheck columns.
                    attrs(dual_profile)["column_units"] = [profile_output_units[data_cols[i]] for i in eachindex(data_cols)]
                    # Place column names
                    attrs(dual_profile)["column_names"] = [data_cols[i] for i in eachindex(data_cols)]
                    # loop over all columns and all zones to store data
                    for i in eachindex(data_cols), k = 1:(sm.props.nz)
                        #storing the partial derivatives with respect to the (partial_index)th input parameter
                        colname = data_cols[i]
                        if colname == "zone" #zone is not a dual number, just copy it
                            dual_profile[k, i] = profile_output_functions[colname](sm, k)
                        else
                            #it happens here: we store the partials
                            #@show colname
                            #@show profile_output_functions[colname](sm, k)
                            #@show profile_output_functions[colname](sm, k).partials
                            #@show profile_output_functions[colname](sm, k).partials[partial_index]
                            dual_profile[k, i] = profile_output_functions[colname](sm, k).partials[partial_index]
                        end
                    end
                    push!(dual_profiles, profile)
                end
                if (!sm.opt.io.hdf5_profile_keep_open)
                    close(sm.profiles_file)
                end
            end
        end
    end
end

function write_terminal_info(sm::StellarModel; now::Bool=false)
    if sm.props.model_number == 1 || sm.props.model_number % sm.opt.io.terminal_header_interval == 0 || now
        print(header)
    end
    if sm.props.model_number == 1 || sm.props.model_number % sm.opt.io.terminal_info_interval == 0 || now
        Printf.format(stdout, line1fmt,
                      sm.props.model_number,
                      log10(sm.props.dt / SECYEAR),
                      log10(get_value(sm.props.L[sm.props.nz])),
                      log10_e * get_value(sm.props.lnT[sm.props.nz]),
                      log10(get_value(sm.props.eos_res[sm.props.nz].P)),
                      log10_e * get_value(sm.props.lnρ[sm.props.nz]),
                      get_value(sm.props.xa[1, sm.network.xa_index[:H1]]),
                      sm.solver_data.newton_iters)
        Printf.format(stdout, line2fmt,
                      sm.props.mstar / MSUN,
                      sm.props.time / SECYEAR,
                      log10_e * get_value(sm.props.lnr[sm.props.nz]) - log10(RSUN),
                      log10_e * get_value(sm.props.lnT[1]),
                      log10(get_value(sm.props.eos_res[1].P)), 
                      log10_e * get_value(sm.props.lnρ[1]),
                      get_value(sm.props.xa[1, sm.network.xa_index[:He4]]),
                      sm.props.nz)
        println()
    end
end


"""
    get_history_dataframe_from_hdf5(hdf5_filename)

Returns a DataFrame object built from an hdf5 file, named `hdf5_filename`.
"""
function get_history_dataframe_from_hdf5(hdf5_filename)
    h5open(hdf5_filename) do history_file
        return DataFrame(history_file["history"][:, :], attrs(history_file["history"])["column_names"])
    end
end

"""
    get_profile_names_from_hdf5(hdf5_filename)

Retruns the column names of the profile data contained in the hdf5 file `hdf5_filename`.
"""
function get_profile_names_from_hdf5(hdf5_filename)
    h5open(hdf5_filename) do profiles_file
        return keys(profiles_file)
    end
end

"""
    get_profile_dataframe_from_hdf5(hdf5_filename, profile_name)

Returns a DataFrame object built from an hdf5 file, named `hdf5_filename`, considering the column named `profile_name`
"""
function get_profile_dataframe_from_hdf5(hdf5_filename, profile_name)
    h5open(hdf5_filename) do profiles_file
        return DataFrame(profiles_file[profile_name][:, :], attrs(profiles_file[profile_name])["column_names"])
    end
end
