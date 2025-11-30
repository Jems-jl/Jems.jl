using ForwardDiff

"""
    mutable struct OneZone{TNUMBER<:Real,TDUALFULL<:ForwardDiff.Dual,
                              TPROPS<:StellarModels.AbstractModelProperties,
                              TNET<:NuclearNetworks.AbstractNuclearNetwork,
                              TSOLVER<:StellarModels.AbstractSolverData}

Structure definition of a model having one internal zone.
"""
@kwdef mutable struct OneZone{TNUMBER<:Real,TDUALFULL<:ForwardDiff.Dual,
                              TPROPS<:StellarModels.AbstractModelProperties,
                              TNET<:NuclearNetworks.AbstractNuclearNetwork,
                              TSOLVER<:StellarModels.AbstractSolverData} <: AbstractModel
    nvars::Int  # This is the sum of hydro vars and species
    var_names::Vector{Symbol}  # List of variable names
    vari::Dict{Symbol,Int}  # Maps variable names to ind_vars vector

    ## Properties related to the solver ##
    # List of equations to be solved.
    composition_equation::StellarModels.TypeStableCompositionEquation{OneZone{TNUMBER,TDUALFULL,TPROPS,TNET,TSOLVER},
                                                                      TDUALFULL}
    # cache to store residuals and solver matrices
    solver_data::TSOLVER

    # Microphyical models
    network::TNET

    # Properties that define the model
    prv_step_props::TPROPS  # properties of the previous step
    props::TPROPS  # properties during and after newton solving

    # Space for used defined options, defaults are in Options.jl
    opt::StellarModels.Options

    # Output files
    history_file::HDF5.File
end

"""
    OneZone(varnames::Vector{Symbol}, composition_equation::Function,
                nvars::Int, nspecies::Int)

Constructor for a `OneZone` instance, using `varnames` for the independent variables, the composition equation
to be solved, number of independent variables `nvars`, number of species in the network `nspecies`
"""
function OneZone(compostion_equation::Function, network::NuclearNetwork, use_static_arrays=true, number_type=Float64)
    nvars = network.nspecies
    var_names_full = network.species_names
    # link var_names to the correct index so you can do ind_var[vari[:lnT]] = 'some temperature'
    vari::Dict{Symbol,Int} = Dict()
    for i in eachindex(network.species_names)
        vari[var_names_full[i]] = i
    end

    solver_data = StellarModels.SolverData(nvars, 1, 0, use_static_arrays, number_type)

    # properties
    prv_step_props = OneZoneProperties(nvars, length(network.reactions), network.nspecies, number_type)
    props = OneZoneProperties(nvars, length(network.reactions), network.nspecies, number_type)

    # create type stable function objects
    dual_sample = ForwardDiff.Dual(zero(number_type), (zeros(number_type, 3 * nvars)...))
    tpe_stbl_func = StellarModels.TypeStableCompositionEquation{OneZone{number_type,typeof(dual_sample),typeof(props),
                                                                        typeof(network),typeof(solver_data)},
                                                                typeof(dual_sample)}(compostion_equation)

    opt = StellarModels.Options()  # create options object

    # create the stellar model
    oz = OneZone(; nvars=nvars,
                 var_names=var_names_full, vari=vari,
                 solver_data=solver_data,
                 composition_equation=tpe_stbl_func, network=network,
                 prv_step_props=prv_step_props, props=props,
                 opt=opt,
                 history_file=HDF5.File(-1, ""))
    init_IO(oz)
    return oz
end

@kwdef mutable struct OneZoneProperties{TN,TDual,TCellDualData} <: AbstractModelProperties
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
    xa::Vector{TCellDualData}   # dim-less
    xa_dual::Vector{TDual}      # only the cell duals wrt itself

    # rates
    rates::Vector{TCellDualData}  # g^-1 s^-1
    rates_dual::Vector{TDual}     # only cell duals wrt itself

    ϵ_nuc::TN
end

function OneZoneProperties(nvars::Int, nrates::Int, nspecies::Int, ::Type{TN}) where {TN<:Real}
    # define the types
    CDDTYPE = CellDualData{nvars + 1,3 * nvars + 1,TN}  # full dual arrays
    TD = typeof(ForwardDiff.Dual(zero(TN), (zeros(TN, nvars))...))  # only the cell duals

    # create the vector containing the independent variables
    ind_vars = zeros(TN, nvars)

    xa_dual = zeros(TD, nvars)
    xa = Vector{CDDTYPE}(undef, nspecies)
    for j = 1:nspecies
        xa[j] = CellDualData(nvars, TN; is_ind_var=true, ind_var_i=nvars - nspecies + j)
    end
    rates_dual = zeros(TD, nrates)
    rates = Vector{CDDTYPE}(undef, nrates)
    for k = 1:nrates
        rates[k] = CellDualData(nvars, TN)
    end

    return OneZoneProperties(; ind_vars=ind_vars, model_number=zero(Int), dt=zero(TN), dt_next=zero(TN), time=zero(TN),
                             xa=xa, xa_dual=xa_dual, rates=rates, rates_dual=rates_dual, T=zero(TN), ρ=zero(TN),
                             ϵ_nuc=zero(TN))
end

"""
    function evaluate_stellar_model_properties!(oz, props::StellarModelProperties{TDual, TCellDualData}) where
        {TDual <: ForwardDiff.Dual, TCellDualData}

Evaluates the stellar model properties `props` from the `ind_vars` array. The goal is to save the 'state' of the
StellarModel so we can easily get properties like rates, eos, opacity values, and retrace if a retry is called.
This does _not_ update the mesh/ind_vars arrays.
"""
function evaluate_one_zone_properties!(oz, props::OneZoneProperties{TN,TDual}) where {TN<:Real,TDual<:ForwardDiff.Dual}
    # update independent variables
    for j = 1:(oz.network.nspecies)
        update_cell_dual_data_value!(props.xa[j], props.ind_vars[oz.nvars - oz.network.nspecies + j])
        props.xa_dual[j] = get_cell_dual(props.xa[j])
    end

    # evaluate rates
    set_rates_for_network!(props.rates_dual, oz.network, props.T, props.ρ, props.xa_dual)
    for j in eachindex(props.rates_dual)
        update_cell_dual_data!(props.rates[j], props.rates_dual[j])
    end

    props.ϵ_nuc = 0.0
    for j in eachindex(props.rates_dual)
        props.ϵ_nuc += props.rates_dual[j].value * oz.network.reactions[j].Qvalue
    end
end

function cycle_props!(oz::OneZone)
    temp_props = oz.prv_step_props
    oz.prv_step_props = oz.props
    oz.props = temp_props
end

function uncycle_props!(oz::OneZone)
    temp_props = oz.props
    oz.props = oz.prv_step_props
    oz.prv_step_props = temp_props
end
