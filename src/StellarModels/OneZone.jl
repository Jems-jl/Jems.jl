using ForwardDiff
using LaTeXStrings

"""
    mutable struct OneZone{TNUMBER<:Real,TDUALFULL<:ForwardDiff.Dual,
                              TPROPS<:StellarModels.AbstractModelProperties,
                              TNET<:NuclearNetworks.AbstractNuclearNetwork,
                              TSOLVER<:StellarModels.AbstractSolverData}

Structure definition of a model having one internal zone.
"""
@kwdef mutable struct OneZone{TPROPS<:StellarModels.AbstractModelProperties,
                              TNET<:NuclearNetworks.AbstractNuclearNetwork,
                              TSOLVER<:StellarModels.AbstractSolverData,
                              TEQUATIONSET<:AbstractEquationSet} <: AbstractModel
    nvars::Int  # This is the sum of hydro vars and species
    var_names::Vector{Symbol}  # List of variable names
    vari::Dict{Symbol,Int}  # Maps variable names to ind_vars vector

    equation_set::TEQUATIONSET

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

    # Output options
    history_output_units::Dict{String,String} = Dict()
    history_output_functions::Dict{String,Function} = Dict()
    history_output_labels::Dict{String,LaTeXStrings.LaTeXString} = Dict()
end

"""
    OneZone(varnames::Vector{Symbol}, composition_equation::Function,
                nvars::Int, nspecies::Int)

Constructor for a `OneZone` instance, using `varnames` for the independent variables, the composition equation
to be solved, number of independent variables `nvars`, number of species in the network `nspecies`
"""
function OneZone(equation_set::AbstractEquationSet, network::NuclearNetwork, use_static_arrays=true,
                    number_type=Float64, internal_dual_tag=ForwardDiff.Tag{:internal, nothing})
    nvars = network.nspecies
    var_names_full = network.species_names
    # link var_names to the correct index so you can do ind_var[vari[:lnT]] = 'some temperature'
    vari::Dict{Symbol,Int} = Dict()
    for i in eachindex(network.species_names)
        vari[var_names_full[i]] = i
    end

    solver_data = build_solver_data_for_equation_set(equation_set, nvars, 1, 0, use_static_arrays, number_type, internal_dual_tag)

    # properties
    prv_step_props = OneZoneProperties(nvars, length(network.reactions), network.nspecies, number_type, internal_dual_tag)
    props = OneZoneProperties(nvars, length(network.reactions), network.nspecies, number_type, internal_dual_tag)

    opt = StellarModels.Options()  # create options object

    # create the stellar model
    oz = OneZone(; nvars=nvars,
                 var_names=var_names_full, vari=vari,
                 equation_set=equation_set,
                 solver_data=solver_data,
                 network=network,
                 prv_step_props=prv_step_props, props=props,
                 opt=opt,
                 history_file=HDF5.File(-1, ""),
                 history_output_units = Dict{String,String}(),
                 history_output_functions = Dict{String,Function}(),
                 history_output_labels = Dict{String,LaTeXStrings.LaTeXString}(),
                 )
    init_IO(oz)
    return oz
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
