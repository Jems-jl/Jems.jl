export DefaultStellarEquationSet, DefaultOneZoneEquationSet

abstract type AbstractEquationSet end

struct DefaultStellarEquationSet<:AbstractEquationSet
end

function hydro_vars(equation_set::DefaultStellarEquationSet)
    return [:lnρ, :lnT, :lnr, :lum]
end

function hydro_vars_scaling(equation_set::DefaultStellarEquationSet)
    return [:log, :log, :log, :maxval]
end

function remesh_splitting(equation_set::DefaultStellarEquationSet, sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
    StellarModels.split_lnr_lnρ(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
    StellarModels.split_lum(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
    StellarModels.split_lnT(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
    StellarModels.split_xa(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
end


function build_properties_for_equation_set(equation_set::DefaultStellarEquationSet, nvars, nz, nextra, network, vari, number_type)
    StellarModelProperties(nvars, nz, nextra,
                                   length(network.reactions), network.nspecies, vari, number_type)
end

struct DefaultOneZoneEquationSet<:AbstractEquationSet
end
