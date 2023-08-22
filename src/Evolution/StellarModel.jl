using BlockBandedMatrices
using SparseArrays
using LinearSolve
using PreallocationTools
using ForwardDiff

"""
    mutable struct StellarStepInfo

Information used for a simulation step. A single stellar model can have three different
objects of type StellarStepInfo, containing information from the previous step, information
right before the Newton solver, and information after the Newton solver has completed.
"""
@kwdef mutable struct StellarStepInfo{T1<:Real}
    # grid properties
    nz::Int # number of zones in the model
    m::Vector{T1} # mass coordinate of each cell
    dm::Vector{T1} # mass contained in each cell
    mstar::T1 # total model mass

    # unique valued properties (ie not cell dependent)
    time::T1
    dt::T1
    model_number::Int

    # full vector with independent variables (size is number of variables * number of 
    # zones)
    ind_vars::Vector{T1}

    # Values of properties at each cell, sizes are equal to number of zones
    lnT::Vector{T1}
    L::Vector{T1}
    lnP::Vector{T1}
    lnρ::Vector{T1}
    lnr::Vector{T1}
end

"""
    mutable struct StellarModel

An evolutionary model for a star, containing information about the star's current
state, as well as the independent variables of the model and its equations.
"""
@kwdef mutable struct StellarModel{T1<:Real, T2<:Real}
    # Properties that define the model
    ind_vars::Vector{T1}  # List of independent variables
    varnames::Vector{Symbol}  # List of variable names
    vari::Dict{Symbol,Int}  # Maps variable names to ind_vars vector
    eqs::Matrix{T2}  # Stores the results of the equation evaluations
    eqs_nums::Vector{T1}
    nvars::Int  # This is the sum of hydro vars and species
    nspecies::Int  # Just the number of species in the network
    structure_equations::Vector{Function}  # List of equations to be solved. Be careful,
    #  typing here probably kills type inference
    diff_caches::Matrix{DiffCache{Vector{T1}, Vector{T1}}}
    # Grid properties
    nz::Int  # Number of zones in the model
    m::Vector{T1}  # Mass coordinate of each cell (g)
    dm::Vector{T1}  # Mass contained in each cell (g)
    mstar::Real  # Total model mass (g)

    # Unique valued properties (ie not cell dependent)
    time::Real  # Age of the model (s)
    dt::Real  # Timestep of the current evolutionary step (s)
    model_number::Int

    # Some basic info
    eos::EOS.AbstractEOS
    eos_results::Matrix{T2}  # size (nz, num_eos_results)
    opacity::Opacity.AbstractOpacity
    opacity_results::Vector{T2}  # size nz
    convection::Convection.AbstractConvection
    conv_results::Matrix{T2}  # size (nz, num_conv_results)
    ∇::Vector{T2}

    isotope_data::Dict{Symbol,Isotope}

    # Jacobian matrix
    jacobian::SparseMatrixCSC{T1, Int64}
    linear_solver::Any  # solver that is produced by LinearSolve

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

    csi::StellarStepInfo  # pointer to current step info appropriate for the phase of the evolution loop 

    # Space for used defined options, defaults are in Options.jl
    opt::Options
end

"""
    StellarModel(varnames::Vector{Symbol}, structure_equations::Vector{Function},
        nvars::Int, nspecies::Int, nz::Int, eos::AbstractEOS, opacity::AbstractOpacity)

Constructor for a `StellarModel` instance, using `varnames` for the independent variables,
functions of the `structure_equations` to be solved, number of independent variables
`nvars`, number of species in the network `nspecies` number of zones in the model
`nz` and an iterface to the EOS and Opacity laws.
"""
function StellarModel(varnames::Vector{Symbol}, structure_equations::Vector{Function},
                      nvars::Int, nspecies::Int, nz::Int,
                      eos::AbstractEOS, opacity::AbstractOpacity, convection::AbstractConvection;
                      solver_method=KLUFactorization())
    dual_type = ForwardDiff.Dual(0.0, (zeros(3 * nvars))...)  # determines type of the AD entries
    ind_vars = zeros(nz * nvars)
    eqs = Matrix{typeof(dual_type)}(undef, nz, nvars)
    eqs_numbers = similar(ind_vars)
    for k = 1:nz
        for i = 1:nvars
            eqs[k, i] = ForwardDiff.Dual(0.0, (zeros(3 * nvars)...))
        end
    end
    dc_type = DiffCache(zeros(nvars), 3 * nvars)
    diff_caches = Matrix{typeof(dc_type)}(undef, nz, 3)
    for k = 1:nz
        for i = 1:3
            diff_caches[k, i] = DiffCache(zeros(nvars), 3 * nvars)
        end
    end
    _init_diff_caches!(diff_caches, nz, nvars)
    opacity_results = Vector{typeof(dual_type)}(undef, nz)
    mlt_results = Matrix{typeof(dual_type)}(undef, nz, convection.num_results)
    eos_results = Matrix{typeof(dual_type)}(undef, nz, eos.num_results)
    ∇ = Vector{typeof(dual_type)}(undef, nz)

    l, u = 1, 1  # block bandwidths
    N = M = nz  # number of row/column blocks
    cols = rows = [nvars for i = 1:N]  # block sizes

    jac_BBM = BlockBandedMatrix(Ones(sum(rows), sum(cols)), rows, cols, (l, u))
    jacobian = sparse(jac_BBM)
    #create solver
    problem = LinearProblem(jacobian, eqs_numbers)
    linear_solver = init(problem, solver_method)

    isotope_data = Chem.get_isotope_list()

    vari::Dict{Symbol,Int} = Dict()
    for i in eachindex(varnames)
        vari[varnames[i]] = i
    end

    dm = zeros(nz)
    m = zeros(nz)

    psi = StellarStepInfo(nz=nz, m=zeros(nz), dm=zeros(nz), mstar=0.0, time=0.0, dt=0.0, model_number=0,
                          ind_vars=ind_vars, lnT=zeros(nz), L=zeros(nz), lnP=zeros(nz),
                          lnρ=zeros(nz), lnr=zeros(nz))
    ssi = StellarStepInfo(nz=nz, m=zeros(nz), dm=zeros(nz), mstar=0.0, time=0.0, dt=0.0, model_number=0,
                          ind_vars=ind_vars, lnT=zeros(nz), L=zeros(nz), lnP=zeros(nz),
                          lnρ=zeros(nz), lnr=zeros(nz))
    esi = StellarStepInfo(nz=nz, m=zeros(nz), dm=zeros(nz), mstar=0.0, time=0.0, dt=0.0, model_number=0,
                          ind_vars=ind_vars, lnT=zeros(nz), L=zeros(nz), lnP=zeros(nz),
                          lnρ=zeros(nz), lnr=zeros(nz))

    opt = Options()

    StellarModel(ind_vars=ind_vars, varnames=varnames, eqs=eqs, eqs_nums=eqs_numbers, nvars=nvars, nspecies=nspecies,
                 structure_equations=structure_equations, diff_caches=diff_caches, vari=vari, nz=nz, m=m, dm=dm,
                 mstar=0.0, time=0.0, dt=0.0,
                 model_number=0, eos=eos, eos_results=eos_results, opacity=opacity,
                 opacity_results=opacity_results, convection=convection, conv_results=mlt_results,
                 ∇=∇, isotope_data=isotope_data, jacobian=jacobian, linear_solver=linear_solver, psi=psi,
                 ssi=ssi, esi=esi, csi=psi, opt=opt)
end

function _init_diff_caches!(dc, nz, nvars)
    for k = 1:nz
        for i = 1:nvars  # set values in `du`
            dc[k, 1].du[i] = 0.0
            dc[k, 2].du[i] = 0.0
            dc[k, 3].du[i] = 0.0
        end
        for i = 1:3  # set all dual values in `dual_du` to 0 for now.
            # dual_du has size (2/3 * nvars * nvars to hold dx_i/dx_j for its own cell + neighbors)
            dc[k, 1].dual_du[:] .= 0.0
            dc[k, 2].dual_du[:] .= 0.0
            dc[k, 3].dual_du[:] .= 0.0
        end
    end
    for k = 1:nz
        for i = 1:nvars
            # these indices are headache inducing
            # diff_caches[k][2].dual_du has structure:
            # (value, dx_1^k/dx_1^k-1, ..., dx_1^k/dx_n^k-1, dx_1^k/dx_1^k, ..., dx_1^k/dx_n^k, dx_1^k/dx_1^k+1, ..., dx_1^k/dx_n^k+1,  # subsize 3*nvars+1
            #  ...
            #  value, dx_n^k/dx_1^k-1, ..., dx_n^k/dx_n^k-1, dx_n^k/dx_1^k, ..., dx_n^k/dx_n^k, dx_n^k/dx_1^k+1, ..., dx_n^k/dx_n^k+1)
            # diff_caches[k][3].dual_du has numerators k -> k+1
            # diff_caches[k][1].dual_du has numerators k -> k-1
            dc[k, 1].dual_du[(i - 1) * (3 * nvars + 1) + 1 + i] = 1.0
            dc[k, 2].dual_du[(i - 1) * (3 * nvars + 1) + 1 + nvars + i] = 1.0
            dc[k, 3].dual_du[(i - 1) * (3 * nvars + 1) + 1 + 2 * nvars + i] = 1.0
        end
    end
end
