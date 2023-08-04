# Perhaps later add support for HDF5, this would come in the form of a function that transforms a history.data into HDF5
# using HDF5
using Printf

function history_get_ind_vars_edge_value(sm::StellarModel, var_symbol::Symbol, edge::Symbol)
    if var_symbol ∉ sm.varnames
        throw(ArgumentError(":$var_symbol is not a valid independent variable"))
    end 
    if edge==:center
        return sm.ind_vars[sm.vari[var_symbol]]
    elseif edge==:surface
        return sm.ind_vars[(sm.nz-1)*sm.nvars + sm.vari[var_symbol]]
    else
        throw(ArgumentError("'edge' must be either :surface or :center. Instead received :$edge"))
    end
end

history_output_options = Dict(
    #general properties
    :star_age => ("year", "star_age", sm->sm.esi.time/SECYEAR),
    :dt => ("year", "dt", sm->sm.esi.dt/SECYEAR),
    :model_number => ("unitless", "model_number", sm->sm.esi.model_number),
    :star_mass => ("Msun", "star_mass", sm->sm.esi.mstar/MSUN),

    #surface properties
    :R_surf => ("Rsun", "R_surf", sm->exp(sm.esi.lnr[sm.nz])/RSUN),
    :L_surf => ("Lsun", "L_surf", sm->sm.esi.L[sm.nz]),
    :T_surf => ("K", "T_surf", sm->exp(sm.esi.lnT[sm.nz])),
    :P_surf => ("dyne", "P_surf", sm->exp(sm.esi.lnP[sm.nz])),
    :ρ_surf => ("g*cm^-3", "ρ_surf", sm->exp(sm.esi.lnρ[sm.nz])),
    :X_surf => ("unitless", "X_surf", sm->history_get_ind_vars_edge_value(sm, :H1, :surface)),
    :Y_surf => ("unitless", "Y_surf", sm->history_get_ind_vars_edge_value(sm, :He4, :surface)),

    #central properties
    :T_center => ("K", "T_center", sm->exp(sm.esi.lnT[1])),
    :P_center => ("dyne", "P_center", sm->exp(sm.esi.lnP[1])),
    :ρ_center => ("g*cm^-3", "ρ_center", sm->exp(sm.ssi.lnρ[1])),
    :X_center => ("unitless", "X_center", sm->history_get_ind_vars_edge_value(sm, :H1, :surface)),
    :Y_center => ("unitless", "Y_center", sm->history_get_ind_vars_edge_value(sm, :He4, :surface)),
)

function write_history_data(sm)
    if (sm.model_number% sm.opt.io.history_interval == 0)
        if !isdir(sm.opt.io.LOGS_directory)
            mkdir(sm.opt.io.LOGS_directory)
        end
        file_existed = isfile(sm.opt.io.LOGS_directory*"/history.data")
        open(sm.opt.io.LOGS_directory*"/history.data","a") do history_file
            data_cols = sm.opt.io.history_values
            if !file_existed
                # need to create the header
                # simply star by creating a row with numbers for each column
                str = ""
                format_string = "%$(sm.opt.io.column_size).0i"
                for i in eachindex(data_cols)
                    if data_cols[i] ∉ keys(StellarEvolution.history_output_options)
                        throw(ArgumentError("Invalid name for history data column, :$(data_cols[i])"))
                    end
                    str = str*Printf.format(Printf.Format(format_string), i)
                end
                println(history_file, str)

                # next up, include the units for all quantities. No need to recheck columns.
                str = ""
                format_string = "%$(sm.opt.io.column_size)s"
                for i in eachindex(data_cols)
                    str = str*Printf.format(Printf.Format(format_string), StellarEvolution.history_output_options[data_cols[i]][1])
                end
                println(history_file, str)

                # Finally, place column names
                str = ""
                for i in eachindex(data_cols)
                    str = str*Printf.format(Printf.Format(format_string), StellarEvolution.history_output_options[data_cols[i]][2])
                end
                println(history_file, str)
            end

            # after being sure the header is there, print the data
            str = ""
            format_string = "%$(sm.opt.io.column_size).$(sm.opt.io.decimal_places)e"
            for i in eachindex(data_cols)
                str = str*Printf.format(Printf.Format(format_string), StellarEvolution.history_output_options[data_cols[i]][3](sm))
            end
            println(history_file, str)
        end
    end
end