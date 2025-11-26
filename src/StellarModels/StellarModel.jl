using ForwardDiff
using LinearAlgebra
using HDF5
using Jems.DualSupport

"""
    mutable struct StellarModel{TN<:Real,TD<:Real,TEOS<:EOS.AbstractEOS,TKAP<:Opacity.AbstractOpacity}

An evolutionary model for a star, containing information about the star's current state, as well as the independent
variables of the model and its equations.

The struct has four parametric types, `TN` for 'normal' numbers, `TD` for dual numbers used in automatic
differentiation, `TEOS` for the type of EOS being used and `TKAP` for the type of opacity law being used.
"""
@kwdef mutable struct StellarModel{TNUMBER<:Real,TDUALFULL<:ForwardDiff.Dual,TPROPS<:AbstractModelProperties,
                                   TEOS<:EOS.AbstractEOS,TKAP<:Opacity.AbstractOpacity,
                                   TNET<:NuclearNetworks.AbstractNuclearNetwork,TTURB<:Turbulence.AbstractTurb,
                                   TSOLVER<:AbstractSolverData} <: AbstractModel
    nvars::Int  # This is the sum of hydro vars and species
    var_names::Vector{Symbol}  # List of variable names
    var_scaling::Vector{Symbol}
    vari::Dict{Symbol,Int}  # Maps variable names to ind_vars vector
    nextra::Int  # Number of extra zones used to avoid constant reallocation while remeshing

    ## Properties related to the solver ##
    # original vector of functions that are solved. These are turned into TypeStableEquations.
    # We keep the original input for when we resize the stellar model.
    structure_equations_original::Vector{Function}
    composition_equation_original::Function
    # List of equations to be solved.
    structure_equations::Vector{TypeStableStructureEquation{StellarModel{TNUMBER,TDUALFULL,TPROPS,TEOS,TKAP,TNET,TTURB,
                                                                         TSOLVER},TDUALFULL}}
    composition_equation::TypeStableCompositionEquation{StellarModel{TNUMBER,TDUALFULL,TPROPS,TEOS,TKAP,TNET,TTURB,
                                                                     TSOLVER},TDUALFULL}
    # cache to store residuals and solver matrices
    solver_data::TSOLVER

    # Remeshing functions
    remesh_split_functions::Vector{Function}

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
end

"""
    StellarModel(varnames::Vector{Symbol}, structure_equations::Vector{Function},
                nvars::Int, nspecies::Int, nz::Int, eos::AbstractEOS, opacity::AbstractOpacity)

Constructor for a `StellarModel` instance, using `varnames` for the independent variables, functions of the
`structure_equations` to be solved, number of independent variables `nvars`, number of species in the network `nspecies`
number of zones in the model `nz` and an iterface to the EOS and Opacity laws.
"""
function StellarModel(var_names::Vector{Symbol}, var_scaling::Vector{Symbol},
                      structure_equations::Vector{Function}, composition_equation::Function, nz::Int, nextra::Int,
                      remesh_split_functions::Vector{Function},
                      network::NuclearNetwork, eos::AbstractEOS, opacity::AbstractOpacity, turbulence::AbstractTurb;
                      use_static_arrays=true, number_type=Float64)
    nvars = length(var_names) + network.nspecies

    # var_names should also contain the name of species, we get them from the network
    var_names_full = vcat(var_names, network.species_names)
    var_scaling_full = vcat(var_scaling, [:unity for i in 1:network.nspecies])

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
    dual_sample = ForwardDiff.Dual(zero(number_type), (zeros(number_type, 3 * nvars)...))
    tpe_stbl_funcs = Vector{TypeStableStructureEquation{StellarModel{number_type,typeof(dual_sample),typeof(props),
                                                                     typeof(eos),typeof(opacity),typeof(network),
                                                                     typeof(turbulence),typeof(solver_data)},
                                                        typeof(dual_sample)}}(undef, length(structure_equations))
    for i in eachindex(structure_equations)
        tpe_stbl_funcs[i] = TypeStableStructureEquation{StellarModel{number_type,typeof(dual_sample),typeof(props),
                                                                     typeof(eos),typeof(opacity),typeof(network),
                                                                     typeof(turbulence),typeof(solver_data)},
                                                        typeof(dual_sample)}(structure_equations[i])
    end
    stbl_comp_eq = TypeStableCompositionEquation{StellarModel{number_type,typeof(dual_sample),typeof(props),
                                                              typeof(eos),typeof(opacity),typeof(network),
                                                              typeof(turbulence),typeof(solver_data)},
                                                              typeof(dual_sample)}(composition_equation)

    opt = Options()  # create options object

    # create the stellar model
    sm = StellarModel(;nvars=nvars,
                      var_names=var_names_full, var_scaling=var_scaling_full,
                      vari=vari, nextra=nextra,
                      solver_data = solver_data,
                      structure_equations_original=structure_equations,
                      structure_equations=tpe_stbl_funcs,
                      composition_equation_original=composition_equation,
                      composition_equation=stbl_comp_eq,
                      remesh_split_functions=remesh_split_functions,
                      eos=eos, opacity=opacity, network=network, turbulence=turbulence,
                      start_step_props=start_step_props, prv_step_props=prv_step_props, props=props,
                      opt=opt,
                      history_file=HDF5.File(-1, ""),
                      profiles_file=HDF5.File(-1, ""))
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

"""
    adjusted_stellar_model_data(sm, new_nz::Int, new_nextra::Int)

Returns a new copy of sm with an adjusted allocated size. This creates a full duplicate
without removing the old stellar model, which is not very memory friendly. One
possible optimization for the future. The new model is created to have `new_nz`
zones with an extra padding of `new_nextra` zones to allow for remeshing.
The new model will copy the contents of

  - ind_vars
  - mstar
  - m
  - dm
  - time
  - dt
  - model_number
  - psi, ssi, esi
  - opt
    As well as the nuclear network, opacity and EOS.
"""
function adjusted_stellar_model_data(sm, new_nz::Int, new_nextra::Int)
    # verify that new size can contain old sm
    if sm.nz > new_nz + new_nextra
        throw(ArgumentError("Can't fit model of size nz=$(sm.nz) using new_nz=$(new_nz) and new_nextra=$(new_nextra)."))
    end
    #get var_names without species
    var_names = sm.var_names[1:(sm.nvars - sm.network.nspecies)]

    new_sm = StellarModel(var_names, sm.structure_equations_original, sm.composition_equation_original,
                          new_nz, new_nextra, sm.remesh_split_functions,
                          sm.network, sm.eos, sm.opacity, sm.turbulence)
    new_sm.nz = sm.nz # If this needs to be adjusted it will be done by remeshing routines
    new_sm.opt = sm.opt

    # backup scalar quantities
    new_sm.time = sm.time
    new_sm.dt = sm.dt
    new_sm.model_number = sm.model_number
    new_sm.mstar = sm.mstar

    new_sm.history_file = sm.history_file
    new_sm.profiles_file = sm.profiles_file

    # copy arrays
    for i = 1:(sm.nz)
        for j = 1:(sm.nvars)
            new_sm.ind_vars[(i - 1) * sm.nvars + j] = sm.ind_vars[(i - 1) * sm.nvars + j]
        end
        new_sm.m[i] = sm.m[i]
        new_sm.dm[i] = sm.dm[i]
    end

    # Copy StellarStepInfo objects
    for (new_ssi, old_ssi) in [(new_sm.psi, sm.psi), (new_sm.ssi, sm.ssi), (new_sm.esi, sm.esi)]
        new_ssi.nz = old_ssi.nz
        new_ssi.time = old_ssi.time
        new_ssi.dt = old_ssi.dt
        new_ssi.model_number = old_ssi.model_number
        new_ssi.mstar = old_ssi.mstar
        for i = 1:(sm.nz)
            for j = 1:(sm.nvars)
                new_ssi.ind_vars[(i - 1) * sm.nvars + j] = old_ssi.ind_vars[(i - 1) * sm.nvars + j]
            end
            new_ssi.m[i] = old_ssi.m[i]
            new_ssi.dm[i] = old_ssi.dm[i]
            new_ssi.lnT[i] = old_ssi.lnT[i]
            new_ssi.L[i] = old_ssi.L[i]
            new_ssi.lnP[i] = old_ssi.lnP[i]
            new_ssi.lnρ[i] = old_ssi.lnρ[i]
            new_ssi.lnr[i] = old_ssi.lnr[i]
            new_ssi.X[i] == old_ssi.X[i]
            new_ssi.Y[i] == old_ssi.Y[i]
        end
    end

    return new_sm
end
