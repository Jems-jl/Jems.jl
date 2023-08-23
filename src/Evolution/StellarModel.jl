using BlockBandedMatrices
using SparseArrays
using LinearSolve
using PreallocationTools

"""
    mutable struct StellarStepInfo

Information used for a simulation step. A single stellar model can have three different objects of type StellarStepInfo,
containing information from the previous step, information right before the Newton solver, and information after
the Newton solver has completed.

The struct has one parametric type, T1 to represent 'normal' numbers. No fields here need to have dual numbers as these
will not be used in automatic differentiation routines.
"""
@kwdef mutable struct StellarStepInfo{T1<:Real}
    # grid properties
    nz::Int  # number of zones in the model
    m::Vector{T1}  # mass coordinate of each cell
    dm::Vector{T1}  # mass contained in each cell
    mstar::T1  # total model mass

    # unique valued properties (ie not cell dependent)
    time::T1
    dt::T1
    model_number::Int

    # full vector with independent variables (size is number of variables * number of zones)
    ind_vars::Vector{T1}

    # Values of properties at each cell, sizes are equal to number of zones
    lnT::Vector{T1}
    L::Vector{T1}
    lnP::Vector{T1}
    lnρ::Vector{T1}
    lnr::Vector{T1}
end

"""
    mutable struct StellarModel{T1<:Real, T2<: Real}

An evolutionary model for a star, containing information about the star's current state, as well as the independent
variables of the model and its equations.

The struct has two parametric types, `T1` for 'normal' numbers, `T2` for dual numbers used in automatic differentiation
"""
@kwdef mutable struct StellarModel{T1<:Real,T2<:Real}
    # Properties that define the model
    ind_vars::Vector{T1}  # List of independent variables
    nvars::Int  # This is the sum of hydro vars and species
    var_names::Vector{Symbol}  # List of variable namesv
    vari::Dict{Symbol,Int}  # Maps variable names to ind_vars vector

    nspecies::Int  # Just the number of species in the network
    species_names::Vector{Symbol}  # just the species names

    # Properties related to the solver
    structure_equations::Vector{Function}  # List of equations to be solved. Be careful,
    #  typing here probably kills type inference
    eqs_numbers::Vector{T1}  # Stores the results of the equation evaluations (as numbers), size nz * nvars
    eqs_duals::Matrix{T2}  # Stores the dual results of the equation evaluation, shape (nz, nvars)
    diff_caches::Matrix{DiffCache{Vector{T1},Vector{T1}}}  # Allocates space for when automatic differentiation needs
    # to happen
    jacobian::SparseMatrixCSC{T1,Int64}  # Jacobian matrix
    linear_solver::Any  # solver that is produced by LinearSolve

    # Grid properties
    nz::Int  # Number of zones in the model
    m::Vector{T1}  # Mass coordinate of each cell (g)
    dm::Vector{T1}  # Mass contained in each cell (g)
    mstar::T1  # Total model mass (g)

    # Unique valued properties (ie not cell dependent)
    time::T1  # Age of the model (s)
    dt::T1  # Timestep of the current evolutionary step (s)
    model_number::Int

    # Some basic info
    eos::EOS.AbstractEOS
    opacity::Opacity.AbstractOpacity

    # Here I want to preemt things that will be necessary once we have an adaptative
    # mesh. Idea is that psi contains the information from the previous step (or the
    # initial condition). ssi will contain information after remeshing. Absent remeshing
    # it will make no difference. esi will contain properties once the step is completed.
    # Information coming from the previous step (psi=Previous Step Info)
    psi::StellarStepInfo
    # Information computed at the start of the step (ssi=Start Step Info)
    ssi::StellarStepInfo
    # Information computed at the end of the step (esi=End Step Info)
    esi::StellarStepInfo

    # Space for used defined options, defaults are in Options.jl
    opt::Options
end

"""
    StellarModel(varnames::Vector{Symbol}, structure_equations::Vector{Function},
        nvars::Int, nspecies::Int, nz::Int, eos::AbstractEOS, opacity::AbstractOpacity)

Constructor for a `StellarModel` instance, using `varnames` for the independent variables, functions of the
`structure_equations` to be solved, number of independent variables `nvars`, number of species in the network `nspecies`
number of zones in the model `nz` and an iterface to the EOS and Opacity laws.
"""
function StellarModel(var_names::Vector{Symbol}, structure_equations::Vector{Function}, nvars::Int, nspecies::Int,
                      nz::Int, eos::AbstractEOS, opacity::AbstractOpacity; solver_method=KLUFactorization())
    ind_vars = ones(nvars * nz)
    species_names = var_names[(nvars - nspecies + 1):end]

    eqs_numbers = ones(nvars * nz)

    dual_sample = ForwardDiff.Dual(0.0, (zeros(3 * nvars)...))
    eqs_duals = Matrix{typeof(dual_sample)}(undef, nz, nvars)
    for k = 1:nz
        for i = 1:nvars
            eqs_duals[k, i] = ForwardDiff.Dual(0.0, (zeros(3 * nvars)...))
        end
    end

    dc_type = DiffCache(zeros(nvars), 3 * nvars)
    diff_caches = Matrix{typeof(dc_type)}(undef, nz, 3)
    for k = 1:nz
        for i = 1:3
            diff_caches[k, i] = DiffCache(zeros(nvars), 3 * nvars)
        end
    end
    l, u = 1, 1  # block bandwidths
    N = M = nz  # number of row/column blocks
    cols = rows = [nvars for i = 1:N]  # block sizes

    jac_BBM = BlockBandedMatrix(Ones(sum(rows), sum(cols)), rows, cols, (l, u))
    jacobian = sparse(jac_BBM)
    #  create solver
    problem = LinearProblem(jacobian, eqs_numbers)
    linear_solver = init(problem, solver_method)

    vari::Dict{Symbol,Int} = Dict()
    for i in eachindex(var_names)
        vari[var_names[i]] = i
    end

    dm = zeros(nz)
    m = zeros(nz)

    psi = StellarStepInfo(nz=nz, m=zeros(nz), dm=zeros(nz), mstar=0.0, time=0.0, dt=0.0, model_number=0,
                          ind_vars=zeros(nvars * nz), lnT=zeros(nz), L=zeros(nz), lnP=zeros(nz), lnρ=zeros(nz),
                          lnr=zeros(nz))
    ssi = StellarStepInfo(nz=nz, m=zeros(nz), dm=zeros(nz), mstar=0.0, time=0.0, dt=0.0, model_number=0,
                          ind_vars=zeros(nvars * nz), lnT=zeros(nz), L=zeros(nz), lnP=zeros(nz), lnρ=zeros(nz),
                          lnr=zeros(nz))
    esi = StellarStepInfo(nz=nz, m=zeros(nz), dm=zeros(nz), mstar=0.0, time=0.0, dt=0.0, model_number=0,
                          ind_vars=zeros(nvars * nz), lnT=zeros(nz), L=zeros(nz), lnP=zeros(nz), lnρ=zeros(nz),
                          lnr=zeros(nz))

    opt = Options()

    StellarModel(ind_vars=ind_vars, var_names=var_names, species_names=species_names, eqs_numbers=eqs_numbers,
                 eqs_duals=eqs_duals, nvars=nvars, nspecies=nspecies, structure_equations=structure_equations,
                 diff_caches=diff_caches, vari=vari, nz=nz, m=m, dm=dm, mstar=0.0, time=0.0, dt=0.0, model_number=0,
                 eos=eos, opacity=opacity, jacobian=jacobian, linear_solver=linear_solver, psi=psi, ssi=ssi, esi=esi,
                 opt=opt)
end
