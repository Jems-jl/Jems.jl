using PreallocationTools
using FunctionWrappers
using StaticArrays

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
struct TypeStableEquation{TS,TD<:Real}
    func::FunctionWrappers.FunctionWrapper{TD,
                                           Tuple{TS,Int,
                                                 Matrix{TD},Matrix{TD},Matrix{TD},
                                                 EOSResults{TD},EOSResults{TD},EOSResults{TD},
                                                 Matrix{TD},
                                                 TD,TD,TD}}
end

"""
    mutable struct StellarStepInfo{TN<:Real}

Information used for a simulation step. A single stellar model can have three different objects of type StellarStepInfo,
containing information from the previous step, information right before the Newton solver, and information after
the Newton solver has completed.

The struct has one parametric type, TN to represent 'normal' numbers. No fields here need to have dual numbers as these
will not be used in automatic differentiation routines.
"""
@kwdef mutable struct StellarStepInfo{TN<:Real}
    # grid properties
    nz::Int  # number of zones in the model
    m::Vector{TN}  # mass coordinate of each cell
    dm::Vector{TN}  # mass contained in each cell
    mstar::TN  # total model mass

    # unique valued properties (ie not cell dependent)
    time::TN
    dt::TN
    model_number::Int

    # full vector with independent variables (size is number of variables * number of zones)
    ind_vars::Vector{TN}

    # Values of properties at each cell, sizes are equal to number of zones
    lnT::Vector{TN}
    L::Vector{TN}
    lnP::Vector{TN}
    lnρ::Vector{TN}
    lnr::Vector{TN}

    #eos results
    eos_res::Vector{EOSResults{TN}}
end

"""
    mutable struct StellarModel{TN<:Real,TD<:Real,TEOS<:EOS.AbstractEOS,TKAP<:Opacity.AbstractOpacity}

An evolutionary model for a star, containing information about the star's current state, as well as the independent
variables of the model and its equations.

The struct has four parametric types, `TN` for 'normal' numbers, `TD` for dual numbers used in automatic
differentiation, `TEOS` for the type of EOS being used and `TKAP` for the type of opacity law being used.
"""
@kwdef mutable struct StellarModel{TN<:Real,TD<:Real,
                                   TEOS<:EOS.AbstractEOS,TKAP<:Opacity.AbstractOpacity,
                                   TSM<:AbstractMatrix,TSV<:AbstractVector}
    # Properties that define the model
    ind_vars::Vector{TN}  # List of independent variables
    nvars::Int  # This is the sum of hydro vars and species
    var_names::Vector{Symbol}  # List of variable namesv
    vari::Dict{Symbol,Int}  # Maps variable names to ind_vars vector

    # Properties related to the solver
    structure_equations_original::Vector{Function} # original vector of functions that are solved. These are turned into TypeStableEquations. We keep the original input for when we resize the stellar model.
    structure_equations::Vector{TypeStableEquation{StellarModel{TN,TD,TEOS,TKAP},TD}}  # List of equations to be solved.
    eqs_numbers::Vector{TN}  # Stores the results of the equation evaluations (as numbers), size nz * nvars
    eqs_duals::Matrix{TD}  # Stores the dual results of the equation evaluation, shape (nz, nvars)
    diff_caches::Matrix{DiffCache{Vector{TN},Vector{TN}}}  # Allocates space for when automatic differentiation needs
    # to happen
    jacobian_D::Vector{TSM}
    jacobian_U::Vector{TSM}
    jacobian_L::Vector{TSM}
    jacobian_tmp::Vector{TSM}
    solver_β::Vector{TSV}
    solver_x::Vector{TSV}
    solver_corr::Vector{TN}

    # Grid properties
    nz::Int  # Number of zones in the model
    nextra::Int  # Number of extra zones used to avoid constant reallocation while remeshing
    m::Vector{TN}  # Mass coordinate of each cell (g)
    dm::Vector{TN}  # Mass contained in each cell (g)
    mstar::TN  # Total model mass (g)

    # Unique valued properties (ie not cell dependent)
    time::TN  # Age of the model (s)
    dt::TN  # Timestep of the current evolutionary step (s)
    model_number::Int

    # Some basic info
    eos::TEOS
    opacity::TKAP
    network::NuclearNetwork

    # cache for the EOS
    eos_res::Matrix{EOSResults{TD}}

    # cache for the rates
    rates_res::Matrix{TD}

    varp1::Matrix{TD}
    var00::Matrix{TD}
    varm1::Matrix{TD}

    # Here I want to preemt things that will be necessary once we have an adaptative
    # mesh. Idea is that psi contains the information from the previous step (or the
    # initial condition). ssi will contain information after remeshing. Absent remeshing
    # it will make no difference. esi will contain properties once the step is completed.
    # Information coming from the previous step (psi=Previous Step Info)
    psi::StellarStepInfo{TN}
    # Information computed at the start of the step (ssi=Start Step Info)
    ssi::StellarStepInfo{TN}
    # Information computed at the end of the step (esi=End Step Info)
    esi::StellarStepInfo{TN}

    # Space for used defined options, defaults are in Options.jl
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
function StellarModel(var_names::Vector{Symbol}, structure_equations::Vector{Function},
                      nz::Int, nextra::Int,
                      network::NuclearNetwork, eos::AbstractEOS, opacity::AbstractOpacity)
    nvars = length(var_names) + network.nspecies

    # create the vector containing the independent variables
    ind_vars = zeros(nvars * (nz + nextra))

    # create the equation results matrix, holding dual numbers (for automatic differentiation, AD)
    dual_sample = ForwardDiff.Dual(0.0, (zeros(3 * nvars)...))
    eqs_duals = Matrix{typeof(dual_sample)}(undef, nz+nextra, nvars)
    for k = 1:(nz + nextra)
        for i = 1:nvars
            eqs_duals[k, i] = ForwardDiff.Dual(0.0, (zeros(3 * nvars)...))
        end
    end

    # create the diff caches
    dc_type = DiffCache(zeros(nvars), 3 * nvars)
    diff_caches = Matrix{typeof(dc_type)}(undef, nz+nextra, 3)
    for k = 1:(nz+nextra)
        for i = 1:3
            diff_caches[k, i] = DiffCache(zeros(nvars), 3 * nvars)
        end
    end

    # create jacobian matrix (we have the diagonal and the upper and lower blocks)
    # we use static arrays, provided by StaticArrays. These are faster than regular
    # arrays for small nvars
    jacobian_D = [(@MMatrix zeros(nvars,nvars)) for i=1:(nz+nextra)]
    jacobian_U = [(@MMatrix zeros(nvars,nvars)) for i=1:(nz+nextra)]
    jacobian_L = [(@MMatrix zeros(nvars,nvars)) for i=1:(nz+nextra)]
    jacobian_tmp = [(@MMatrix zeros(nvars,nvars)) for i=1:(nz+nextra)]
    solver_β = [(@MVector zeros(nvars)) for i=1:(nz+nextra)]
    solver_x = [(@MVector zeros(nvars)) for i=1:(nz+nextra)]
    solver_corr = zeros(nvars*(nz+nextra))

    # create the equation results vector for the solver (holds plain numbers instead of duals)
    eqs_numbers = ones(nvars * (nz+nextra))

    # var_names should also contain the name of species, we get them from the network
    var_names_full = vcat(var_names, network.species_names)

    # link var_names to the correct index so you can do ind_var[vari[:lnT]] = 'some temperature'
    vari::Dict{Symbol,Int} = Dict()
    for i in eachindex(var_names_full)
        vari[var_names_full[i]] = i
    end

    # create type stable function objects
    tsfs = Vector{TypeStableEquation{StellarModel{eltype(ind_vars),typeof(dual_sample),typeof(eos),typeof(opacity)},
                                     typeof(dual_sample)}}(undef, length(structure_equations))
    for i in eachindex(structure_equations)
        tsfs[i] = TypeStableEquation{StellarModel{eltype(ind_vars),typeof(dual_sample),typeof(eos),typeof(opacity)},
                                     typeof(dual_sample)}(structure_equations[i])
    end

    # mass coordinates
    dm = zeros(nz+nextra)
    m = zeros(nz+nextra)

    # eos results
    eos_res = [EOSResults{typeof(dual_sample)}() for i = 1:(nz+nextra), j = 1:3]

    # rates results
    rates_res = Matrix{typeof(dual_sample)}(undef, (nz+nextra), length(network.reactions))

    # create stellar step info objects
    psi = StellarStepInfo(nz=nz, m=zeros(nz+nextra), dm=zeros(nz+nextra), mstar=0.0, time=0.0, dt=0.0, model_number=0,
                          ind_vars=zeros(nvars * (nz+nextra)), lnT=zeros(nz+nextra), L=zeros(nz+nextra), lnP=zeros(nz+nextra), lnρ=zeros(nz+nextra),
                          lnr=zeros(nz+nextra), eos_res=[EOSResults{Float64}() for i = 1:(nz+nextra)])
    ssi = StellarStepInfo(nz=nz, m=zeros(nz+nextra), dm=zeros(nz+nextra), mstar=0.0, time=0.0, dt=0.0, model_number=0,
                          ind_vars=zeros(nvars * (nz+nextra)), lnT=zeros(nz+nextra), L=zeros(nz+nextra), lnP=zeros(nz+nextra), lnρ=zeros(nz+nextra),
                          lnr=zeros(nz+nextra), eos_res=[EOSResults{Float64}() for i = 1:(nz+nextra)])
    esi = StellarStepInfo(nz=nz, m=zeros(nz+nextra), dm=zeros(nz+nextra), mstar=0.0, time=0.0, dt=0.0, model_number=0,
                          ind_vars=zeros(nvars * (nz+nextra)), lnT=zeros(nz+nextra), L=zeros(nz+nextra), lnP=zeros(nz+nextra), lnρ=zeros(nz+nextra),
                          lnr=zeros(nz+nextra), eos_res=[EOSResults{Float64}() for i = 1:(nz+nextra)])

    # create options object
    opt = Options()

    plt = Plotter(nothing)


    # create the stellar model
    sm = StellarModel(ind_vars=ind_vars, var_names=var_names_full, eqs_numbers=eqs_numbers,
                      eqs_duals=eqs_duals, nvars=nvars,
                      structure_equations_original=structure_equations,
                      structure_equations=tsfs,
                      diff_caches=diff_caches, vari=vari, nz=nz, nextra=nextra,
                      m=m, dm=dm, mstar=0.0, time=0.0, dt=0.0, model_number=0,
                      varp1=Matrix{typeof(dual_sample)}(undef, nz+nextra, nvars),
                      var00=Matrix{typeof(dual_sample)}(undef, nz+nextra, nvars),
                      varm1=Matrix{typeof(dual_sample)}(undef, nz+nextra, nvars),
                      eos=eos, opacity=opacity, network=network,
                      jacobian_D=jacobian_D, jacobian_U=jacobian_U, jacobian_L=jacobian_L,
                      jacobian_tmp=jacobian_tmp, solver_β=solver_β,
                      solver_x=solver_x, solver_corr=solver_corr,
                      eos_res=eos_res, rates_res = rates_res,
                      psi=psi, ssi=ssi, esi=esi, opt=opt, plt=plt)
    init_diff_cache!(sm)
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
                      new_nz, new_nextra, sm.network, sm.eos, sm.opacity)
    new_sm.nz = sm.nz # If this needs to be adjusted it will be done by remeshing routines
    new_sm.opt = sm.opt

    # backup scalar quantities
    new_sm.time = sm.time
    new_sm.dt = sm.dt
    new_sm.model_number = sm.model_number
    new_sm.mstar = sm.mstar

    # copy arrays
    for i in 1:sm.nz
        for j in 1:sm.nvars
            new_sm.ind_vars[(i-1)*sm.nvars + j] = sm.ind_vars[(i-1)*sm.nvars + j]
        end
        new_sm.m[i] = sm.m[i]
        new_sm.dm[i] = sm.dm[i]
    end

    # Copy StellarStepInfo objects
    for (new_ssi, old_ssi) in [(new_sm.psi, sm.psi),(new_sm.ssi, sm.ssi),(new_sm.esi, sm.esi)]
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
        end
    end

    return new_sm
end

"""
    init_diff_cache!(sm::StellarModel)

Initializes the diff_caches to the values of the independent variables, and sets ones in the correct spots where the
dx_i^k/dx_i^k entries lie.
"""
function init_diff_cache!(sm::StellarModel)
    for k = 1:(sm.nz)
        # set all partials to 0 for the moment
        sm.diff_caches[k, 1].dual_du[:] .= 0.0
        sm.diff_caches[k, 2].dual_du[:] .= 0.0
        sm.diff_caches[k, 3].dual_du[:] .= 0.0
        #= these indices are headache inducing...
        diff_caches[k, 2].dual_du has structure:
        (x_1, dx_1^k/dx_1^k-1, ..., dx_1^k/dx_n^k-1, dx_1^k/dx_1^k, ..., dx_1^k/dx_n^k, dx_1^k/dx_1^k+1, ..., dx_1^k/dx_n^k+1,  # subsize 3*nvars+1
         ...
         x_n, dx_n^k/dx_1^k-1, ..., dx_n^k/dx_n^k-1, dx_n^k/dx_1^k, ..., dx_n^k/dx_n^k, dx_n^k/dx_1^k+1, ..., dx_n^k/dx_n^k+1)
        diff_caches[k, 3].dual_du has numerators k -> k+1
        diff_caches[k, 1].dual_du has numerators k -> k-1
        =#
        for i = 1:(sm.nvars)
            # set variable values in du, dual_du and its corresponding non-zero derivative
            if k != 1
                sm.diff_caches[k, 1].du[i] = sm.ind_vars[sm.nvars * (k - 2) + i]
                sm.diff_caches[k, 1].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1] = sm.ind_vars[sm.nvars * (k - 2) + i]
                sm.diff_caches[k, 1].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1 + i] = 1.0  # dx^k-1_i/dx^k-1_i = 1!!
            end
            sm.diff_caches[k, 2].du[i] = sm.ind_vars[sm.nvars * (k - 1) + i]
            sm.diff_caches[k, 2].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1] = sm.ind_vars[sm.nvars * (k - 1) + i]
            sm.diff_caches[k, 2].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1 + sm.nvars + i] = 1.0  # dx^k_i/dx^k_i = 1!!
            if k != sm.nz
                sm.diff_caches[k, 3].du[i] = sm.ind_vars[sm.nvars * k + i]
                sm.diff_caches[k, 3].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1] = sm.ind_vars[sm.nvars * k + i]
                sm.diff_caches[k, 3].dual_du[(i - 1) * (3 * sm.nvars + 1) + 1 + 2 * sm.nvars + i] = 1.0  # dx^k+1_i/dx^k+1_i = 1!!
            end
        end
    end
end
