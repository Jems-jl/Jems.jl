using ForwardDiff
using LinearAlgebra
using HDF5
using Jems.DualSupport
using LaTeXStrings

"""
    mutable struct StellarModel{TN<:Real,TD<:Real,TEOS<:EOS.AbstractEOS,TKAP<:Opacity.AbstractOpacity}

An evolutionary model for a star, containing information about the star's current state, as well as the independent
variables of the model and its equations.

The struct has four parametric types, `TN` for 'normal' numbers, `TD` for dual numbers used in automatic
differentiation, `TEOS` for the type of EOS being used and `TKAP` for the type of opacity law being used.
"""
@kwdef mutable struct StellarModel{TPROPS<:AbstractModelProperties,
                                   TEOS<:EOS.AbstractEOS,TKAP<:Opacity.AbstractOpacity,
                                   TNET<:NuclearNetworks.AbstractNuclearNetwork,TTURB<:Turbulence.AbstractTurb,
                                   TSOLVER<:AbstractSolverData, TEQUATIONSET<:AbstractEquationSet} <: AbstractModel
    nvars::Int  # This is the sum of hydro vars and species
    var_names::Vector{Symbol}  # List of variable names
    var_scaling::Vector{Symbol}
    vari::Dict{Symbol,Int}  # Maps variable names to ind_vars vector
    nextra::Int  # Number of extra zones used to avoid constant reallocation while remeshing

    equation_set::TEQUATIONSET

    solver_data::TSOLVER

    # Microphyical models
    eos::TEOS
    opacity::TKAP
    network::TNET
    turbulence::TTURB

    # Properties that define the model
    prv_step_props::TPROPS  # properties of the previous step
    start_step_props::TPROPS  # properties before newton solving (but after remesh)
    props::TPROPS  # properties during and after newton solving

    # Space for used defined options, defaults are in Options.jl
    opt::Options

    # Output files
    history_file::HDF5.File
    profiles_file::HDF5.File

    # Output options
    history_output_units::Dict{String,String} = Dict()
    history_output_functions::Dict{String,Function} = Dict()
    history_output_labels::Dict{String,LaTeXStrings.LaTeXString} = Dict()
    profile_output_units::Dict{String,String} = Dict()
    profile_output_functions::Dict{String,Function} = Dict()
    profile_output_labels::Dict{String,Union{LaTeXStrings.LaTeXString, String}} = Dict()
end

"""
    StellarModel(varnames::Vector{Symbol}, structure_equations::Vector{Function},
                nvars::Int, nspecies::Int, nz::Int, eos::AbstractEOS, opacity::AbstractOpacity)

Constructor for a `StellarModel` instance, using `varnames` for the independent variables, functions of the
`structure_equations` to be solved, number of independent variables `nvars`, number of species in the network `nspecies`
number of zones in the model `nz` and an iterface to the EOS and Opacity laws.
"""
function StellarModel(equation_set::AbstractEquationSet,
                      nz::Int, nextra::Int,
                      network::NuclearNetwork, eos::AbstractEOS, opacity::AbstractOpacity, turbulence::AbstractTurb;
                      use_static_arrays=true, number_type=Float64, internal_dual_tag=ForwardDiff.Tag{:internal, nothing})
    hydro_var_names = hydro_vars(equation_set)
    nvars = length(hydro_var_names) + network.nspecies

    # var_names should also contain the name of species, we get them from the network
    var_names_full = vcat(hydro_var_names, network.species_names)
    var_scaling_full = vcat(hydro_vars_scaling(equation_set), [:unity for i in 1:network.nspecies])

    # link var_names to the correct index so you can do ind_var[vari[:lnT]] = 'some temperature'
    vari::Dict{Symbol,Int} = Dict()
    for i in eachindex(var_names_full)
        vari[var_names_full[i]] = i
    end

    solver_data = build_solver_data_for_equation_set(equation_set, nvars, nz, nextra, use_static_arrays, number_type, internal_dual_tag)

    # properties
    prv_step_props = build_properties_for_equation_set(equation_set, nvars, nz, nextra, network, vari, number_type, internal_dual_tag)
    start_step_props = build_properties_for_equation_set(equation_set, nvars, nz, nextra, network, vari, number_type, internal_dual_tag)
    props = build_properties_for_equation_set(equation_set, nvars, nz, nextra, network, vari, number_type, internal_dual_tag)

    opt = Options()  # create options object

    # create the stellar model
    sm = StellarModel(;nvars=nvars,
                      var_names=var_names_full, var_scaling=var_scaling_full,
                      vari=vari, nextra=nextra,
                      equation_set=equation_set,
                      solver_data = solver_data,
                      eos=eos, opacity=opacity, network=network, turbulence=turbulence,
                      start_step_props=start_step_props, prv_step_props=prv_step_props, props=props,
                      opt=opt,
                      history_file=HDF5.File(-1, ""),
                      profiles_file=HDF5.File(-1, ""),
                      history_output_units = Dict{String,String}(),
                      history_output_functions = Dict{String,Function}(),
                      history_output_labels = Dict{String,LaTeXStrings.LaTeXString}(),
                      profile_output_units = Dict{String,String}(),
                      profile_output_functions = Dict{String,Function}(),
                      profile_output_labels = Dict{String,Union{LaTeXStrings.LaTeXString, String}}())
    init_IO(sm) # initializes output options
    return sm
end

"""
    cycle_props!(sm::StellarModel)

Moves the model properties of the StellarModel `sm` over one state:
start_step_props -> props -> prv_step_props -> start_step_props
"""
function cycle_props!(sm::StellarModel)
    temp_props = sm.prv_step_props
    sm.prv_step_props = sm.props
    sm.props = sm.start_step_props
    sm.start_step_props = temp_props
end

"""
    uncycle_props!(sm::StellarModel)

Moves the model properties of the StellarModel `sm` back one state:
start_step_props <- props <- prv_step_props <- start_step_props
"""
function uncycle_props!(sm::StellarModel)
    temp_props = sm.props
    sm.props = sm.prv_step_props
    sm.prv_step_props = sm.start_step_props
    sm.start_step_props = temp_props
end
