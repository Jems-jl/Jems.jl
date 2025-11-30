using HDF5
using DataFrames
using Printf
using LaTeXStrings

export add_history_option, add_profile_option,
       history_output_units, history_output_functions, history_output_labels,
       profile_output_units, profile_output_functions, profile_output_labels

ioinit = false

const width = 9
const decimals = 4
const floatstr = "%#$width.$decimals" * "g "
const intstr = "%$width" * "i "

mutable struct TerminalHeader
    header::String
    linefmts::Vector{Printf.Format}
end

const terminal_header = TerminalHeader("", [])

function setup_header(sm::StellarModel)
    terminal_header.linefmts = Vector{String}(undef, 2)
    terminal_header.linefmts[1] = Printf.Format(intstr * floatstr^6 * intstr * "\n")
    terminal_header.linefmts[2] = Printf.Format(floatstr^7 * intstr * "\n")

    terminal_header.header = """
        model     logdt      logL   logTeff     logPs     logρs    H_cntr     iters
         mass       age      logR     logTc     logPc     logρc   He_cntr     zones
    -------------------------------------------------------------------------------
    """
end

function setup_header(oz::OneZone)
    lines = oz.network.nspecies ÷ 6 + 1
    lastline = oz.network.nspecies % 6
    if lastline == 0
        terminal_header.linefmts = Vector{Printf.Format}(undef, lines)
    else
        terminal_header.linefmts = Vector{Printf.Format}(undef, lines + 1)
    end
    terminal_header.linefmts[1] = Printf.Format(intstr * floatstr^4 * intstr * "\n")

    j = oz.network.nspecies
    i = 2
    while j > 6
        terminal_header.linefmts[i] = Printf.Format(floatstr^6 * "\n")
        j -= 6
        i += 1
    end
    terminal_header.linefmts[end] = Printf.Format(floatstr^j * "\n")

    terminal_header.header = """
        model     logdt       age      logT      logρ     iters
    """

    for (j, species) in enumerate(oz.network.species_names)
        if j % 6 != 1
            speciesstr = lpad(String(species), width+1)
        else
            speciesstr = lpad(String(species), width)
        end
        terminal_header.header *= speciesstr
        if j == oz.network.nspecies || j % 6 == 0
            terminal_header.header *= "\n"
        end
    end

    terminal_header.header *= """
    -----------------------------------------------------------
    """
end

const history_output_units::Dict{String,String} = Dict()
const history_output_functions::Dict{String,Function} = Dict()
const history_output_labels::Dict{String,LaTeXStrings.LaTeXString} = Dict()
const profile_output_units::Dict{String,String} = Dict()
const profile_output_functions::Dict{String,Function} = Dict()
const profile_output_labels::Dict{String,Union{LaTeXStrings.LaTeXString, String}} = Dict()

function add_history_option(name, unit, func; label::Union{LaTeXStrings.LaTeXString, String}="")
    if haskey(history_output_units, name)
        throw(ArgumentError("Key $name is already part of the history output options"))
    end
    history_output_units[name] = unit
    history_output_functions[name] = func
    if label == ""
        history_output_labels[name] = name
    else
        history_output_labels[name] = label
    end
end

function setup_model_history_functions(sm::StellarModel)
    # general properties
    add_history_option("age", "year", sm -> sm.props.time / SECYEAR, label=L"\text{age}\,[\text{yr}]")
    add_history_option("dt", "year", sm -> sm.props.dt / SECYEAR, label=L"\Delta t\,[\text{yr}]")
    add_history_option("model_number", "unitless", sm -> sm.props.model_number, label=L"\text{Model Number}")
    add_history_option("star_mass", "Msun", sm -> sm.props.mstar / MSUN, label=L"\text{Mass}\,[M_\odot]")

    # surface properties
    add_history_option("R_surf", "Rsun", sm -> exp(get_value(sm.props.lnr[sm.props.nz])) / RSUN, label=L"R_\text{surf}\,[R_\odot]")
    add_history_option("L_surf", "Lsun", sm -> get_value(sm.props.L[sm.props.nz]), label=L"\text{surf}\,[L_\odot]")
    add_history_option("T_surf", "K", sm -> exp(get_value(sm.props.lnT[sm.props.nz])), label=L"T_\text{surf}\,[\text{K}]")
    add_history_option("rho_surf", "g*cm^-3", sm -> exp(get_value(sm.props.lnρ[sm.props.nz])), label=L"\rho_\text{surf}\,[\text{g\,cm^{-3}}]")
    add_history_option("P_surf", "dyne", sm -> exp(get_value(sm.props.eos_res[sm.props.nz].P)), label=L"P_\text{surf}\,[\text{dyne}]")
    add_history_option("X_surf", "unitless", sm -> get_value(sm.props.xa[sm.props.nz, sm.network.xa_index[:H1]]), label=L"X_\text{surf}")
    add_history_option("Y_surf", "unitless", sm -> get_value(sm.props.xa[sm.props.nz, sm.network.xa_index[:He4]]), label=L"Y_\text{surf}")

    # central properties
    add_history_option("T_center", "K", sm -> exp(get_value(sm.props.lnT[1])), label=L"T_\text{c}\,[\text{K}]")
    add_history_option("rho_center", "g*cm^-3", sm -> exp(get_value(sm.props.lnρ[1])), label=L"\rho_\text{c}\,[\text{g\,cm^{-3}}]")
    add_history_option("P_center", "dyne", sm -> get_value(sm.props.eos_res[1].P), label=L"P_\text{c}\,[\text{dyne}]")
    add_history_option("X_center", "unitless", sm -> get_value(sm.props.xa[1, sm.network.xa_index[:H1]]), label=L"X_\text{c}")
    add_history_option("Y_center", "unitless", sm -> get_value(sm.props.xa[1, sm.network.xa_index[:He4]]), label=L"Y_\text{c}")
end

function setup_model_history_functions(oz::OneZone)
    # general properties
    add_history_option("age", "year", oz -> oz.props.time / SECYEAR, label=L"\text{Age}\,[\text{yr}]")
    add_history_option("dt", "year", oz -> oz.props.dt / SECYEAR, label=L"\Delta t\,[\text{yr}]")
    add_history_option("model_number", "unitless", oz -> oz.props.model_number, label="\text{Model Number}")

    add_history_option("T", "K", oz -> oz.props.T, label=L"T\,[\text{K}]")
    add_history_option("rho", "g*cm^-3", oz -> oz.props.ρ, label=L"\rho\,[\text{g\,cm^{-3}}]")
    for j in eachindex(oz.network.species_names)
        species = oz.network.species_names[j]
        add_history_option(String(species), "unitless", oz -> get_value(oz.props.xa[oz.network.xa_index[species]]))
    end
end

function add_profile_option(name, unit, func; label::Union{LaTeXStrings.LaTeXString, String}="")
    if haskey(profile_output_units, name)
        throw(ArgumentError("Key $name is already part of the history output options"))
    end
    profile_output_units[name] = unit
    profile_output_functions[name] = func
    if label == ""
        profile_output_labels[name] = name
    else
        profile_output_labels[name] = label
    end
end

function setup_model_profile_functions(sm::StellarModel)
    # general properties
    add_profile_option("zone", "unitless", (sm, k) -> k, label=L"\text{Zone}")
    add_profile_option("mass", "Msun", (sm, k) -> sm.props.m[k] / MSUN, label=L"\text{Mass}\,[M_\odot]")
    add_profile_option("dm", "Msun", (sm, k) -> sm.props.dm[k] / MSUN, label=L"\Delta m\,[M_\odot]")

    # thermodynamic properties
    add_profile_option("log10_r", "log10(Rsun)", (sm, k) -> get_value(sm.props.lnr[k]) * log10(ℯ) - log10(RSUN), label=L"\log_{10}(r/R_\odot)")
    add_profile_option("log10_P", "log10(dyne)", (sm, k) -> log10(get_value(sm.props.eos_res[k].P)), label=L"\log_{10}(P/[\text{dyne}])")
    add_profile_option("log10_T", "log10(K)", (sm, k) -> get_value(sm.props.lnT[k]) * log10(ℯ), label=L"\log_{10}(T/[\text{K}])")
    add_profile_option("log10_rho", "log10_(g*cm^-3)", (sm, k) -> get_value(sm.props.lnρ[k]) * log10(ℯ), label=L"\log_{10}(\rho/[\text{g\,cm^{-3}}])")
    add_profile_option("luminosity", "Lsun", (sm, k) -> get_value(sm.props.L[k]) / LSUN, label=L"L/L_\odot")

    # abundances
    add_profile_option("X", "unitless", (sm, k) -> get_value(sm.props.xa[k, sm.network.xa_index[:H1]]), label=L"X")
    add_profile_option("Y", "unitless", (sm, k) -> get_value(sm.props.xa[k, sm.network.xa_index[:He4]]), label=L"Y")

    # temperature gradients
    add_profile_option("nabla_a_face", "unitless", (sm, k) -> get_value(sm.props.∇ₐ_face[k]), label=L"\nabla_\text{a,face}")
    add_profile_option("nabla_r_face", "unitless", (sm, k) -> get_value(sm.props.∇ᵣ_face[k]), label=L"\nabla_\text{r,face}")
    add_profile_option("nabla_face", "unitless", (sm, k) -> get_value(sm.props.turb_res[k].∇), label=L"\nabla_\text{face}")
    add_profile_option("D_face", "cm^2*s^{-1}", (sm, k) -> get_value(sm.props.turb_res[k].D_turb), label=L"D_\text{face}\,[\text{cm^2\,s^{-1}}]")
end

function init_IO(m::AbstractModel)
    setup_header(m)
    setup_model_history_functions(m)
    if isa(m, OneZone)
        global ioinit = true
        return
    end
    setup_model_profile_functions(m)
    global ioinit = true
end

function clear_IO()
    for i in eachindex(history_output_functions)
        delete!(history_output_functions, i)
        delete!(history_output_units, i)
    end
    for i in eachindex(profile_output_functions)
        delete!(profile_output_functions, i)
        delete!(profile_output_units, i)
    end
    terminal_header.header = ""
    terminal_header.linefmts = []
    global ioinit = false
end

"""
    create_output_files(sm::StellarModel)

Creates output files for history and profile data
"""
function create_output_files!(m::AbstractModel)
    if !ioinit
        init_IO(m)
    end

    # Create history file
    m.history_file = h5open(m.opt.io.hdf5_history_filename, "w")
    data_cols = m.opt.io.history_values
    ncols = length(data_cols)

    # verify validity of column names
    for i in eachindex(data_cols)
        if data_cols[i] ∉ keys(history_output_functions)
            throw(ArgumentError("Invalid name for history data column, :$(data_cols[i])"))
        end
    end

    # Create history dataset in HDF5 file
    # Dataset is created with size (0, ncols), we will add rows by using the HDF5.set_extent_dims function
    # the (-1, ncols) is used to define the maximum extent of the dataset, -1 indicates that it is unbound
    # in number of rows. The chunk size is used for compression. Smaller chunk sizes will result in worse
    # compression but faster writes.
    # The compression level can be anywhere between 0 and 9, 0 being no compression 9 being the highest.
    # Compression is lossless.
    history = create_dataset(m.history_file, "history", Float64, ((0, ncols), (-1, ncols)),
                             chunk=(m.opt.io.hdf5_history_chunk_size, ncols),
                             compress=m.opt.io.hdf5_history_compression_level)

    # next up, include the units for all quantities. No need to recheck columns.
    attrs(history)["column_units"] = [history_output_units[data_cols[i]] for i in eachindex(data_cols)]
    # Finally, place column names
    attrs(history)["column_names"] = [data_cols[i] for i in eachindex(data_cols)]
    if (!m.opt.io.hdf5_history_keep_open)
        close(m.history_file)
    end

    if isa(m, OneZone)
        return
    end

    # Create profile file
    m.profiles_file = h5open(m.opt.io.hdf5_profile_filename, "w")
    data_cols = m.opt.io.profile_values
    # verify validity of column names
    for i in eachindex(data_cols)
        if data_cols[i] ∉ keys(profile_output_functions)
            throw(ArgumentError("Invalid name for profile data column, :$(data_cols[i])"))
        end
    end
    if (!m.opt.io.hdf5_profile_keep_open)
        close(m.profiles_file)
    end
end

function shut_down_IO!(m)
    if (m.opt.io.hdf5_history_keep_open)
        close(m.history_file)
    end
    if (m.opt.io.hdf5_profile_keep_open)
        close(m.profiles_file)
    end
    if ioinit
        clear_IO()
    end
end

"""
    write_data(sm::StellarModel)

Saves data (history/profile) for the current model, as required by the settings in `sm.opt.io`.
"""
function write_data(m::AbstractModel)
    # do history
    if (m.opt.io.history_interval > 0)
        file_exists = isfile(m.opt.io.hdf5_history_filename)
        if !file_exists
            throw(ErrorException("History file does not exist at $(m.opt.io.hdf5_history_filename)"))
        end
        if (m.props.model_number % m.opt.io.history_interval == 0)
            if (!m.opt.io.hdf5_history_keep_open)
                m.history_file = h5open(m.opt.io.hdf5_history_filename, "r+")
            end
            data_cols = m.opt.io.history_values
            ncols = length(data_cols)

            # after being sure the header is there, print the data
            history = m.history_file["history"]
            HDF5.set_extent_dims(history, (size(history)[1] + 1, ncols))
            for i in eachindex(data_cols)
                history[end, i] = history_output_functions[data_cols[i]](m)
            end
            if (!m.opt.io.hdf5_history_keep_open)
                close(m.history_file)
            end
        end
    end

    if isa(m, OneZone)
        return
    end

    # do profile
    if (m.opt.io.profile_interval > 0)
        file_exists = isfile(m.opt.io.hdf5_profile_filename)
        if !file_exists  # create file if it doesn't exist yet
            throw(ErrorException("Profile file does not exist at $(m.opt.io.hdf5_profile_filename)"))
        end
        if (m.props.model_number % m.opt.io.profile_interval == 0)
            if (!m.opt.io.hdf5_profile_keep_open)
                m.profiles_file = h5open(m.opt.io.hdf5_profile_filename, "r+")
            end
            data_cols = m.opt.io.profile_values
            ncols = length(data_cols)
            # Save current profile
            profile = create_dataset(m.profiles_file,
                                     "$(lpad(m.props.model_number,m.opt.io.hdf5_profile_dataset_name_zero_padding,"0"))",
                                     Float64, ((m.props.nz, ncols), (m.props.nz, ncols));
                                     chunk=(m.opt.io.hdf5_profile_chunk_size, ncols),
                                     compress=m.opt.io.hdf5_profile_compression_level)

            # next up, include the units for all quantities. No need to recheck columns.
            attrs(profile)["column_units"] = [profile_output_units[data_cols[i]] for i in eachindex(data_cols)]
            # Place column names
            attrs(profile)["column_names"] = [data_cols[i] for i in eachindex(data_cols)]

            # store data
            for i in eachindex(data_cols), k = 1:(m.props.nz)
                profile[k, i] = profile_output_functions[data_cols[i]](m, k)
            end
            if (!m.opt.io.hdf5_profile_keep_open)
                close(m.profiles_file)
            end
        end
    end
end

function write_terminal_info(sm::StellarModel; now::Bool=false)
    if sm.props.model_number == 1 || sm.props.model_number % sm.opt.io.terminal_header_interval == 0 || now
        print(terminal_header.header)
    end
    if sm.props.model_number == 1 || sm.props.model_number % sm.opt.io.terminal_info_interval == 0 || now
        Printf.format(stdout, terminal_header.linefmts[1],
                      sm.props.model_number,
                      log10(sm.props.dt / SECYEAR),
                      log10(get_value(sm.props.L[sm.props.nz])),
                      log10(ℯ) * get_value(sm.props.lnT[sm.props.nz]),
                      log10(get_value(sm.props.eos_res[sm.props.nz].P)),
                      log10(ℯ) * get_value(sm.props.lnρ[sm.props.nz]),
                      get_value(sm.props.xa[1, sm.network.xa_index[:H1]]),
                      sm.solver_data.newton_iters)
        Printf.format(stdout, terminal_header.linefmts[2],
                      sm.props.mstar / MSUN,
                      sm.props.time / SECYEAR,
                      log10(ℯ) * get_value(sm.props.lnr[sm.props.nz]) - log10(RSUN),
                      log10(ℯ) * get_value(sm.props.lnT[1]),
                      log10(get_value(sm.props.eos_res[1].P)),
                      log10(ℯ) * get_value(sm.props.lnρ[1]),
                      get_value(sm.props.xa[1, sm.network.xa_index[:He4]]),
                      sm.props.nz)
        println()
    end
end

function write_terminal_info(oz::OneZone; now::Bool=false)
    if oz.props.model_number == 1 || oz.props.model_number % oz.opt.io.terminal_header_interval == 0 || now
        print(terminal_header.header)
    end
    if oz.props.model_number == 1 || oz.props.model_number % oz.opt.io.terminal_info_interval == 0 || now
        Printf.format(stdout, terminal_header.linefmts[1],
                      oz.props.model_number,
                      log10(oz.props.dt / SECYEAR),
                      oz.props.time / SECYEAR,
                      log10(oz.props.T),
                      log10(oz.props.ρ),
                      oz.solver_data.newton_iters)
        j = oz.network.nspecies
        i = 1
        while j > 6
            Printf.format(stdout, terminal_header.linefmts[i + 1], (get_value.(oz.props.xa[((i - 1) * 6 + 1):(6i)]))...)
            i += 1
            j -= 6
        end
        Printf.format(stdout, terminal_header.linefmts[end], (get_value.(oz.props.xa[((i - 1) * 6 + 1):end]))...)
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
