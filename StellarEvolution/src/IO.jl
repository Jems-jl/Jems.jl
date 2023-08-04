# Perhaps later add support for HDF5, this would come in the form of a function that transforms a history.data into HDF5
# using HDF5

function history_get_ind_vars_edge_value(sm::StellarModel, var_symbol::Symbol, edge::Symbol)
    if var_symbol ∉ sm.vari
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
    :star_mass => ("Msun", "star_mass", sm->sm.esi.mass/MSUN),

    #surface properties
    :R_surf => ("Rsun", "R_surf", sm->exp(sm.esi.lnr[sm.nz])/RSUN),
    :L_surf => ("Lsun", "L_surf", sm->sm.esi.L[sm.nz]/LSUN),
    :T_surf => ("K", "T_surf", sm->exp(sm.esi.lnT[sm.nz])),
    :P_surf => ("dyne", "P_surf", sm->exp(sm.esi.lnP[sm.nz])),
    :ρ_surf => ("g*cm^-3", "ρ_surf", sm->exp(sm.esi.lnρ[sm.nz])),
    :X_surf => ("unitless", "X_surf", sm->exp(history_get_ind_vars_edge_value(sm, :H1, :surface))),
    :Y_surf => ("unitless", "Y_surf", sm->exp(history_get_ind_vars_edge_value(sm, :He4, :surface))),

    #central properties
    :T_center => ("K", "T_center", sm->exp(sm.esi.lnT[1])),
    :P_center => ("dyne", "P_center", sm->exp(sm.esi.lnP[1])),
    :ρ_center => ("g*cm^-3", "ρ_center", sm->exp(sm.ssi.lnρ[1])),
    :X_center => ("unitless", "X_center", sm->exp(history_get_ind_vars_edge_value(sm, :lnT, :surface))),
    :Y_center => ("unitless", "Y_center", sm->exp(history_get_ind_vars_edge_value(sm, :He4, :surface))),
)

