using FunctionWrappers
using ForwardDiff
using LinearAlgebra
using HDF5
using Jems.DualSupport

"""
    struct TypeStableEquation{TS,TD<:Real}

Structure that wraps a stellar structure equation into a type stable object, using FunctionWrappers.jl. This requires
that the stellar structure equations have the following signature:

    ```
    function structure_equation(::TS, ::Int)::TD
    ```

For typical usage, TS is the concrete type of StellarModel, and TD the type of dual number being used for automatic
differentiation. The function must return an object of type TD, the result of the equation.
"""
struct TypeStableEquation{TS,TD<:Real}
    func::FunctionWrappers.FunctionWrapper{TD,Tuple{TS,Int}}
end

"""

    mutable struct StellarModel{TN<:Real,TD<:Real,TEOS<:EOS.AbstractEOS,TKAP<:Opacity.AbstractOpacity}

An evolutionary model for a star, containing information about the star's current state, as well as the independent
variables of the model and its equations.

The struct has many parametric types, `TN` for 'normal' numbers, `TD` for dual numbers of `TN`, used in automatic
differentiation, `TPROPS` for the internal properties of the model, `TEOS` for the type of EOS, `TKAP` for the type of
opacity law, `TNET` for the nuclear network, and `TSOLVER` for the structures that solve the linear system.
"""
@kwdef mutable struct StellarModel{TNUMBER<:Real,TDUALFULL<:ForwardDiff.Dual,TPROPS<:AbstractStellarModelProperties,
                                   TEOS<:EOS.AbstractEOS,TKAP<:Opacity.AbstractOpacity,
                                   TNET<:NuclearNetworks.AbstractNuclearNetwork,TSOLVER<:AbstractSolverData}
    # Basic info that does not change over the run (for now)
    nvars::Int  # This is the sum of hydro vars and species
    var_names::Vector{Symbol}  # List of variable names
    vari::Dict{Symbol,Int}  # Maps variable names to ind_vars vector
    nextra::Int  # Number of extra zones used to avoid constant reallocation while remeshing

    ## Properties related to the solver ##
    # original vector of functions that are solved. These are turned into TypeStableEquations.
    # We keep the original input for when we resize the stellar model.
    structure_equations_original::Vector{Function}
    # List of equations to be solved.
    structure_equations::Vector{
        TypeStableEquation{StellarModel{TNUMBER,TDUALFULL,TPROPS,TEOS,TKAP,TNET,TSOLVER},TDUALFULL}}
    # cache to store residuals and solver matrices
    solver_data::TSOLVER

    # Remeshing functions
    remesh_split_functions::Vector{Function}

    # Microphyical models
    eos::TEOS
    opacity::TKAP
    network::TNET

    # Properties that define the model
    prv_step_props::TPROPS  # properties of the previous step
    start_step_props::TPROPS  # properties before newton solving (but after remesh)
    props::TPROPS  # properties during and after newton solving

    # Space for used defined options, defaults are in Options.jl
    opt::Options

    # object holding plotting things, ie figures, data to plot.
    plt::Plotter

    # Output files
    history_file::HDF5.File
    profiles_file::HDF5.File
end

"""

    StellarModel(varnames::Vector{Symbol}, structure_equations::Vector{Function},
                nvars::Int, nspecies::Int, nz::Int, eos::AbstractEOS, opacity::AbstractOpacity)

Constructor for a `StellarModel` instance, using `varnames` for the independent variables, functions of the
`structure_equations` to be solved, number of independent variables `nvars`, number of species in the network `nspecies`
number of zones in the model `nz` and an iterface to the EOS and Opacity laws.
"""
function StellarModel(var_names::Vector{Symbol},
                      structure_equations::Vector{Function}, nz::Int, nextra::Int,
                      remesh_split_functions::Vector{Function},
                      network::NuclearNetwork, eos::AbstractEOS, opacity::AbstractOpacity;
                      use_static_arrays=true, number_type=Float64)
    nvars = length(var_names) + network.nspecies

    # var_names should also contain the name of species, we get them from the network
    var_names_full = vcat(var_names, network.species_names)

    # link var_names to the correct index so you can do ind_var[vari[:lnT]] = 'some temperature'
    vari::Dict{Symbol,Int} = Dict()
    for i in eachindex(var_names_full)
        vari[var_names_full[i]] = i
    end

    solver_data = SolverData(nvars, nz, nextra, use_static_arrays, number_type)

    # properties
    prv_step_props = StellarModelProperties(nvars, nz, nextra,
                                       length(network.reactions), network.nspecies, vari, number_type)
    start_step_props = StellarModelProperties(nvars, nz, nextra,
                                       length(network.reactions), network.nspecies, vari, number_type)
    props = StellarModelProperties(nvars, nz, nextra,
                                       length(network.reactions), network.nspecies, vari, number_type)

    # create type stable function objects
    dual_sample = ForwardDiff.Dual(zero(number_type), (zeros(number_type, 3*nvars)...))
    tpe_stbl_funcs = Vector{TypeStableEquation{StellarModel{number_type,typeof(dual_sample),typeof(props),
                                                            typeof(eos),typeof(opacity),typeof(network),
                                                            typeof(solver_data)},
                                     typeof(dual_sample)}}(undef, length(structure_equations))
    for i in eachindex(structure_equations)
        tpe_stbl_funcs[i] = TypeStableEquation{StellarModel{number_type,typeof(dual_sample),typeof(props),
                                                            typeof(eos),typeof(opacity),typeof(network),
                                                            typeof(solver_data)},
                                     typeof(dual_sample)}(structure_equations[i])
    end

    opt = Options()  # create options object
    plt = Plotter()

    # create the stellar model
    sm = StellarModel(; nvars=nvars,
                      var_names=var_names_full, vari=vari, nextra=nextra,
                      solver_data = solver_data,
                      structure_equations_original=structure_equations,
                      structure_equations=tpe_stbl_funcs,
                      remesh_split_functions=remesh_split_functions,
                      eos=eos, opacity=opacity, network=network,
                      start_step_props=start_step_props, prv_step_props=prv_step_props, props=props,
                      opt=opt, plt=plt,
                      history_file=HDF5.File(-1,""),
                      profiles_file=HDF5.File(-1,""))
    return sm
end

