using FunctionWrappers
using ForwardDiff
using LinearAlgebra
using Jems.DualSupport

"""
    struct TypeStableEquation{TS,TD<:Real}

Structure that wraps a stellar structure equation into a type stable object, using FunctionWrappers.jl. This requires
that the stellar structure equations have the following signature:

```
function structure_equation(::TS, ::Int,
                            ::Matrix{TD}, ::Matrix{TD}, ::Matrix{TD},
                            ::EOSResults{TD}, ::EOSResults{TD}, ::EOSResults{TD}
                            ::Matrix{TD},
                            ::TD, ::TD, ::TD)::TD
```

For typical usage, TS is the concrete type of StellarModel, and TD the type of dual number being used for automatic
differentiation. The function must return an object of type TD, the result of the equation.
"""
struct TypeStableEquation{TPROPS<:AbstractStellarModelProperties,Options,TD<:Real}
    func::FunctionWrappers.FunctionWrapper{TD,
                                           Tuple{TPROPS,Options,Int}}
end

"""
    mutable struct StellarModel{TN<:Real,TD<:Real,TEOS<:EOS.AbstractEOS,TKAP<:Opacity.AbstractOpacity}

An evolutionary model for a star, containing information about the star's current state, as well as the independent
variables of the model and its equations.

The struct has four parametric types, `TN` for 'normal' numbers, `TD` for dual numbers used in automatic
differentiation, `TEOS` for the type of EOS being used and `TKAP` for the type of opacity law being used.
"""
@kwdef mutable struct StellarModel{TDUALFULL<:ForwardDiff.Dual, TPROPS<:AbstractStellarModelProperties,
                                   TEOS<:EOS.AbstractEOS,TKAP<:Opacity.AbstractOpacity,TNET<:NuclearNetworks.AbstractNuclearNetwork,
                                   TTURB<:Turbulence.AbstractTurb, TSOLVER<:AbstractSolverData}
    ## Properties related to the solver ##
    # original vector of functions that are solved. These are turned into TypeStableEquations.
    # We keep the original input for when we resize the stellar model.
    structure_equations_original::Vector{Function}
    # List of equations to be solved.
    structure_equations::Vector{TypeStableEquation{TPROPS,Options,TDUALFULL}}

    solver_data::TSOLVER

    # Remeshing functions
    remesh_split_functions::Vector{Function}

    # Some basic info
    eos::TEOS
    opacity::TKAP
    network::TNET
    turbulence::TTURB

    ##
    props::TPROPS # properties derived from ind_vars
    props_old::TPROPS # properties from previous timestep
    props_old_after_remesh::TPROPS # properties from previous timestep after remeshing

    # Space for user defined options, defaults are in Options.jl
    opt::Options

    # object holding plotting things, ie figures, data to plot.
    plt::Plotter
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
                      network::NuclearNetwork, eos::AbstractEOS, opacity::AbstractOpacity, turbulence::AbstractTurb;
                      use_static_arrays=true, number_type=Float64)
    nvars = length(var_names) + network.nspecies

    # create the vector containing the independent variables

    # var_names should also contain the name of species, we get them from the network
    var_names_full = vcat(var_names, network.species_names)

    solver_data = SolverData(nvars, nz, nextra, use_static_arrays, number_type)

    props = StellarModelProperties(var_names_full, nz, nextra, 
                    length(network.reactions), network.nspecies, number_type)
    props_old = StellarModelProperties(var_names_full, nz, nextra, 
                    length(network.reactions), network.nspecies, number_type)
    props_old_after_remesh = StellarModelProperties(var_names_full, nz, nextra, 
                    length(network.reactions), network.nspecies, number_type)

    # create type stable function objects
    dual_sample = ForwardDiff.Dual(zero(number_type), (zeros(number_type, 3*nvars)...))
    tpe_stbl_funcs = Vector{TypeStableEquation{typeof(props), Options,
                                typeof(dual_sample)}}(undef, length(structure_equations))
    for i in eachindex(structure_equations)
        tpe_stbl_funcs[i] = TypeStableEquation{typeof(props), Options,
                                     typeof(dual_sample)}(structure_equations[i])
    end

    # create options object
    opt = Options()

    plt = Plotter()

    # create the stellar model
    sm = StellarModel(structure_equations_original=structure_equations,
                      structure_equations=tpe_stbl_funcs,
                      solver_data = solver_data,
                      remesh_split_functions=remesh_split_functions,
                      eos=eos, opacity=opacity, network=network, turbulence=turbulence,
                      props=props, props_old=props_old, props_old_after_remesh=props_old_after_remesh,
                      opt=opt, plt=plt)

    return sm
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
    if sm.nz > new_nz+new_nextra
        throw(ArgumentError("Can't fit model of size nz=$(sm.nz) using new_nz=$(new_nz) and new_nextra=$(new_nextra)."))
    end
    #get var_names without species
    var_names = sm.var_names[1:sm.nvars-sm.network.nspecies]

    new_sm = StellarModel(var_names, sm.structure_equations_original,
                      new_nz, new_nextra, sm.remesh_split_functions,
                      sm.network, sm.eos, sm.opacity, sm.turbulence)
    new_sm.nz = sm.nz # If this needs to be adjusted it will be done by remeshing routines
    new_sm.opt = sm.opt

    # backup scalar quantities
    new_sm.time = sm.time
    new_sm.dt = sm.dt
    new_sm.model_number = sm.model_number
    new_sm.mstar = sm.mstar
    new_sm.plt = sm.plt

    # copy arrays
    for i in 1:sm.nz
        for j in 1:sm.nvars
            new_sm.ind_vars[(i-1)*sm.nvars + j] = sm.ind_vars[(i-1)*sm.nvars + j]
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
        for i in 1:sm.nz
            for j in 1:sm.nvars
                new_ssi.ind_vars[(i-1)*sm.nvars + j] = old_ssi.ind_vars[(i-1)*sm.nvars + j]
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
