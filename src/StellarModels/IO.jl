using HDF5
using DataFrames
using Printf

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

"""
    history_get_ind_vars_edge_value(sm::StellarModel, var_symbol::Symbol, edge::Symbol)

Returns the value of the independent variable `var_symbol` at either the surface or the center of the StellarModel `sm`.
`edge` can be either `:surface` or `:center`.
"""
function history_get_ind_vars_edge_value(sm::StellarModel, var_symbol::Symbol, edge::Symbol)
    if var_symbol ∉ sm.var_names
        throw(ArgumentError(":$var_symbol is not a valid independent variable"))
    end
    if edge == :center
        return sm.ind_vars[sm.vari[var_symbol]]
    elseif edge == :surface
        return sm.ind_vars[(sm.nz - 1) * sm.nvars + sm.vari[var_symbol]]
    else
        throw(ArgumentError("'edge' must be either :surface or :center. Instead received 
                :$edge"))
    end
end

history_output_options = Dict(
                              #general properties
                              "star_age" => ("year", sm -> sm.esi.time / SECYEAR),
                              "dt" => ("year", sm -> sm.esi.dt / SECYEAR),
                              "model_number" => ("unitless", sm -> sm.esi.model_number),
                              "star_mass" => ("Msun", sm -> sm.esi.mstar / MSUN),

                              #surface properties
                              "R_surf" => ("Rsun", sm -> exp(sm.esi.lnr[sm.nz]-log(RSUN))),
                              "L_surf" => ("Lsun", sm -> sm.esi.L[sm.nz]),
                              "T_surf" => ("K", sm -> exp(sm.esi.lnT[sm.nz])),
                              "P_surf" => ("dyne", sm -> exp(sm.esi.lnP[sm.nz])),
                              "ρ_surf" => ("g*cm^-3", sm -> exp(sm.esi.lnρ[sm.nz])),
                              "X_surf" => ("unitless", sm -> history_get_ind_vars_edge_value(sm, :H1, :surface)),
                              "Y_surf" => ("unitless", sm -> history_get_ind_vars_edge_value(sm, :He4, :surface)),

                              #central properties
                              "T_center" => ("K", sm -> exp(sm.esi.lnT[1])),
                              "P_center" => ("dyne", sm -> exp(sm.esi.lnP[1])),
                              "ρ_center" => ("g*cm^-3", sm -> exp(sm.ssi.lnρ[1])),
                              "X_center" => ("unitless", sm -> history_get_ind_vars_edge_value(sm, :H1, :surface)),
                              "Y_center" => ("unitless", sm -> history_get_ind_vars_edge_value(sm, :He4, :surface)))

"""
    profile_get_ind_vars_value(sm::StellarModel, var_symbol::Symbol, k::Int)

Returns the value of the variable Symbol `var_symbol` at cell number `k` of the StellarModel `sm`.
"""
function profile_get_ind_vars_value(sm::StellarModel, var_symbol::Symbol, k::Int)
    if var_symbol ∉ sm.var_names
        throw(ArgumentError(":$var_symbol is not a valid independent variable"))
    end
    return sm.ind_vars[(k - 1) * sm.nvars + sm.vari[var_symbol]]
end

profile_output_options = Dict(
                              #general properties
                              "zone" => ("unitless", (sm, k) -> k),
                              "mass" => ("Msun", (sm, k) -> sm.esi.m[k] / MSUN),
                              "dm" => ("Msun", (sm, k) -> sm.esi.dm[k] / MSUN),

                              #thermodynamic properties
                              "log10_r" => ("log10(Rsun)", (sm, k) -> sm.esi.lnr[k] / RSUN * log10_e),
                              "log10_P" => ("log10(dyne)", (sm, k) -> sm.esi.lnP[k] * log10_e),
                              "log10_T" => ("log10(K)", (sm, k) -> sm.esi.lnT[k] * log10_e),
                              "log10_ρ" => ("log10_(g*cm^-3)", (sm, k) -> sm.esi.lnρ[k] * log10_e),
                              "luminosity" => ("Lsun", (sm, k) -> sm.esi.L[k]),

                              #abundance
                              "X" => ("unitless", (sm, k) -> profile_get_ind_vars_value(sm, :H1, k)),
                              "Y" => ("unitless", (sm, k) -> profile_get_ind_vars_value(sm, :He4, k)))

"""
    write_data(dm::StellarModel)

Saves data (history/profile) for the current model, as required by the settings in `sm.opt.io`.
"""
function write_data(sm::StellarModel)
    # do history
    if (sm.opt.io.history_interval > 0)
        file_exists = isfile(sm.opt.io.hdf5_history_filename)
        if !file_exists  # create file if it doesn't exist yet
            h5open(sm.opt.io.hdf5_history_filename, "w") do history_file
                data_cols = sm.opt.io.history_values
                ncols = length(data_cols)

                # verify validity of column names
                for i in eachindex(data_cols)
                    if data_cols[i] ∉ keys(history_output_options)
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
                history = create_dataset(history_file, "history", Float64, ((0, ncols), (-1, ncols)),
                                         chunk=(sm.opt.io.hdf5_history_chunk_size, ncols),
                                         compress=sm.opt.io.hdf5_history_compression_level)

                # next up, include the units for all quantities. No need to recheck columns.
                attrs(history)["column_units"] = [history_output_options[data_cols[i]][1] for i in eachindex(data_cols)]
                # Finally, place column names
                attrs(history)["column_names"] = [data_cols[i] for i in eachindex(data_cols)]
            end
        end
        if (sm.model_number % sm.opt.io.history_interval == 0)
            h5open(sm.opt.io.hdf5_history_filename, "r+") do history_file
                data_cols = sm.opt.io.history_values
                ncols = length(data_cols)

                # after being sure the header is there, print the data
                history = history_file["history"]
                HDF5.set_extent_dims(history, (size(history)[1] + 1, ncols))
                for i in eachindex(data_cols)
                    history[end, i] = history_output_options[data_cols[i]][2](sm)
                end
            end
        end
    end
    # do profile
    if (sm.opt.io.profile_interval > 0)
        file_exists = isfile(sm.opt.io.hdf5_profile_filename)
        if !file_exists  # create file if it doesn't exist yet
            h5open(sm.opt.io.hdf5_profile_filename, "w") do profile_file
                data_cols = sm.opt.io.profile_values
                # verify validity of column names
                for i in eachindex(data_cols)
                    if data_cols[i] ∉ keys(profile_output_options)
                        throw(ArgumentError("Invalid name for history data column,
                            :$(data_cols[i])"))
                    end
                end
            end
        end
        if (sm.model_number % sm.opt.io.profile_interval == 0)
            h5open(sm.opt.io.hdf5_profile_filename, "r+") do profile_file
                data_cols = sm.opt.io.profile_values
                ncols = length(data_cols)
                # Save current profile
                profile = create_dataset(profile_file,
                                         "$(lpad(sm.model_number,sm.opt.io.hdf5_profile_dataset_name_zero_padding,"0"))",
                                         Float64, ((sm.nz, ncols), (sm.nz, ncols));
                                         chunk=(sm.opt.io.hdf5_profile_chunk_size, ncols),
                                         compress=sm.opt.io.hdf5_profile_compression_level)

                # next up, include the units for all quantities. No need to recheck columns.
                attrs(profile)["column_units"] = [profile_output_options[data_cols[i]][1] for i in eachindex(data_cols)]
                # Place column names
                attrs(profile)["column_names"] = [data_cols[i] for i in eachindex(data_cols)]

                # store data
                for i in eachindex(data_cols), k = 1:(sm.nz)
                    profile[k, i] = profile_output_options[data_cols[i]][2](sm, k)
                end
            end
        end
    end
end

function write_terminal_info(sm::StellarModel)
    if sm.esi.model_number == 1 || 
            sm.esi.model_number % sm.opt.io.terminal_header_interval == 0
        print(header)
    end
    if sm.esi.model_number == 1 ||
            sm.esi.model_number % sm.opt.io.terminal_info_interval == 0
        Printf.format(stdout, line1fmt, sm.esi.model_number, log10(sm.dt/SECYEAR), log10(sm.esi.L[sm.nz]),
                      log10_e * sm.esi.lnT[sm.nz], log10_e * sm.esi.lnP[sm.nz], log10_e * sm.esi.lnρ[sm.nz],
                      sm.esi.X[1], sm.solver_data.newton_iters)
        Printf.format(stdout, line2fmt, sm.esi.mstar / MSUN, sm.esi.time / SECYEAR, log10_e * sm.esi.lnr[sm.nz] - log10(RSUN),
                      log10_e * sm.esi.lnT[1], log10_e * sm.esi.lnP[1], log10_e*sm.esi.lnρ[1], sm.esi.Y[1], sm.esi.nz)
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
    h5open(hdf5_filename) do profile_file
        return keys(profile_file)
    end
end

"""
    get_profile_dataframe_from_hdf5(hdf5_filename, profile_name)

Returns a DataFrame object built from an hdf5 file, named `hdf5_filename`, considering the column named `profile_name`
"""
function get_profile_dataframe_from_hdf5(hdf5_filename, profile_name)
    h5open(hdf5_filename) do profile_file
        return DataFrame(profile_file[profile_name][:, :], attrs(profile_file[profile_name])["column_names"])
    end
end
