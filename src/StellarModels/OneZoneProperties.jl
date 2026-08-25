@kwdef mutable struct OneZoneProperties{TN,TDual,TLocalDualData} <: AbstractModelProperties
    # scalar quantities
    dt::TN  # Timestep of the current evolutionary step (s)
    dt_next::TN
    time::TN  # Age of the model (s)
    model_number::Int
    nz = 1  # duh; however  solver needs to know what nz is

    # T and ρ for the zone
    T::TN
    ρ::TN

    # array of the values of the independent variables, everything should be reconstructable from this (mesh dependent)
    ind_vars::Vector{TN}

    # independent variables (duals constructed from the ind_vars array)
    xa::Vector{TLocalDualData}   # dim-less
    xa_dual::Vector{TDual}      # only the cell duals wrt itself

    # rates
    rates::Vector{TLocalDualData}  # g^-1 s^-1
    rates_dual::Vector{TDual}     # only cell duals wrt itself

    ϵ_nuc::TN
end

function OneZoneProperties(nvars::Int, nrates::Int, nspecies::Int,
                            ::Type{TN}, ::Type{internal_tag}) where {TN<:Real, internal_tag<:ForwardDiff.Tag}
    # define the types
    LDDTYPE = LocalDualData{nvars + 1,3 * nvars + 1,TN, internal_tag}  # full dual arrays
    TDL = typeof(ForwardDiff.Dual{internal_tag}(zero(TN), (zeros(TN, nvars))...))  # only the local duals

    # create the vector containing the independent variables
    ind_vars = zeros(TN, nvars)

    xa_dual = zeros(TDL, nvars)
    xa = Vector{LDDTYPE}(undef, nspecies)
    for j = 1:nspecies
        xa[j] = LocalDualData(nvars, TN, internal_tag; is_ind_var=true, ind_var_i=nvars - nspecies + j)
    end
    rates_dual = zeros(TDL, nrates)
    rates = Vector{LDDTYPE}(undef, nrates)
    for k = 1:nrates
        rates[k] = LocalDualData(nvars, TN, internal_tag)
    end

    return OneZoneProperties(; ind_vars=ind_vars, model_number=zero(Int), dt=zero(TN), dt_next=zero(TN), time=zero(TN),
                             xa=xa, xa_dual=xa_dual, rates=rates, rates_dual=rates_dual, T=zero(TN), ρ=zero(TN),
                             ϵ_nuc=zero(TN))
end

"""
    function evaluate_one_zone_properties!(oz, props::OneZoneProperties)

Evaluates the one zone model properties `props` from the `ind_vars` array. The goal is to save the 'state' of the
OneZone so we can easily get properties like rates, eos, opacity values, and retrace if a retry is called.
This does _not_ update the mesh/ind_vars arrays.
"""
function evaluate_one_zone_properties!(oz, props::OneZoneProperties)
    # update independent variables
    for j = 1:(oz.network.nspecies)
        update_local_dual_data_value!(props.xa[j], props.ind_vars[oz.nvars - oz.network.nspecies + j])
        props.xa_dual[j] = get_local_dual(props.xa[j])
    end

    # evaluate rates
    set_rates_for_network!(props.rates_dual, oz.network, props.T, props.ρ, props.xa_dual)
    for j in eachindex(props.rates_dual)
        update_local_dual_data!(props.rates[j], props.rates_dual[j])
    end

    props.ϵ_nuc = 0.0
    for j in eachindex(props.rates_dual)
        props.ϵ_nuc += props.rates_dual[j].value * oz.network.reactions[j].Qvalue
    end
end
