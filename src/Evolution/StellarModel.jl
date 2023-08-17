using BlockBandedMatrices
using SparseArrays
using LinearSolve

"""
    mutable struct StellarStepInfo

Information used for a simulation step. A single stellar model can have three different
objects of type StellarStepInfo, containing information from the previous step, information
right before the Newton solver, and information after the Newton solver has completed.
"""
@kwdef mutable struct StellarStepInfo
    # grid properties
    nz::Int # number of zones in the model
    m::Vector{<:Real} # mass coordinate of each cell
    dm::Vector{<:Real} # mass contained in each cell
    mstar::Real # total model mass

    # unique valued properties (ie not cell dependent)
    time::Real
    dt::Real
    model_number::Int

    # full vector with independent variables (size is number of variables * number of 
    # zones)
    ind_vars::Vector{<:Real}

    # Values of properties at each cell, sizes are equal to number of zones
    lnT::Vector{<:Real}
    L::Vector{<:Real}
    lnP::Vector{<:Real}
    lnρ::Vector{<:Real}
    lnr::Vector{<:Real}
end

"""
    mutable struct StellarModel

An evolutionary model for a star, containing information about the star's current
state, as well as the independent variables of the model and its equations.
"""
@kwdef mutable struct StellarModel
    # Properties that define the model
    ind_vars::Vector{<:Real}  # List of independent variables
    varnames::Vector{Symbol}  # List of variable names
    vari::Dict{Symbol,Int}  # Maps variable names to ind_vars vector
    eqs::Vector{<:Real}  # Stores the results of the equation evaluations
    nvars::Int  # This is the sum of hydro vars and species
    nspecies::Int  # Just the number of species in the network
    structure_equations::Vector{Function}  # List of equations to be solved. Be careful,
    #  typing here probably kills type inference

    # Grid properties
    nz::Int  # Number of zones in the model
    m::Vector{<:Real}  # Mass coordinate of each cell (g)
    dm::Vector{<:Real}  # Mass contained in each cell (g)
    mstar::Real  # Total model mass (g)

    # Unique valued properties (ie not cell dependent)
    time::Real  # Age of the model (s)
    dt::Real  # Timestep of the current evolutionary step (s)
    model_number::Int

    # Some basic info
    eos::EOS.AbstractEOS
    eos_results::Matrix{<:Real}  # size (nz, num_eos_results)
    opacity::Opacity.AbstractOpacity
    opacity_results::Vector{<:Real}  # size nz
    convection::Convection.AbstractConvection
    conv_results::Matrix{<:Real}  # size (nz, num_conv_results)
    ∇::Vector{<:Real}

    isotope_data::Dict{Symbol,Isotope}

    # Jacobian matrix
    jacobian::SparseMatrixCSC{Float64,Int64}
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
    ind_vars = ones(nvars * nz)
    eqs = ones(nvars * nz)
    m = ones(nz)

    l, u = 1, 1  # block bandwidths
    N = M = nz  # number of row/column blocks
    cols = rows = [nvars for i = 1:N]  # block sizes

    jac_BBM = BlockBandedMatrix(Ones(sum(rows), sum(cols)), rows, cols, (l, u))
    jacobian = sparse(jac_BBM)
    #create solver
    problem = LinearProblem(jacobian, eqs)
    linear_solver = init(problem, solver_method)

    isotope_data = Chem.get_isotope_list()

    vari::Dict{Symbol,Int} = Dict()
    for i in eachindex(varnames)
        vari[varnames[i]] = i
    end

    dm = zeros(nz)
    m = zeros(nz)

    psi = StellarStepInfo(nz=nz, m=zeros(nz), dm=zeros(nz), mstar=0.0, time=0.0, dt=0.0, model_number=0,
                          ind_vars=zeros(nvars * nz), lnT=zeros(nz), L=zeros(nz), lnP=zeros(nz),
                          lnρ=zeros(nz), lnr=zeros(nz))
    ssi = StellarStepInfo(nz=nz, m=zeros(nz), dm=zeros(nz), mstar=0.0, time=0.0, dt=0.0, model_number=0,
                          ind_vars=zeros(nvars * nz), lnT=zeros(nz), L=zeros(nz), lnP=zeros(nz),
                          lnρ=zeros(nz), lnr=zeros(nz))
    esi = StellarStepInfo(nz=nz, m=zeros(nz), dm=zeros(nz), mstar=0.0, time=0.0, dt=0.0, model_number=0,
                          ind_vars=zeros(nvars * nz), lnT=zeros(nz), L=zeros(nz), lnP=zeros(nz),
                          lnρ=zeros(nz), lnr=zeros(nz))

    opt = Options()

    StellarModel(ind_vars=ind_vars, varnames=varnames, eqs=eqs, nvars=nvars, nspecies=nspecies,
                 structure_equations=structure_equations, vari=vari, nz=nz, m=m, dm=dm, mstar=0.0, time=0.0, dt=0.0,
                 model_number=0, eos=eos, eos_results=zeros(nz, eos.num_results), opacity=opacity,
                 opacity_results=zeros(nz), convection=convection, conv_results=zeros(nz, convection.num_results), ∇=zeros(nz),
                 isotope_data=isotope_data, jacobian=jacobian, linear_solver=linear_solver, psi=psi, ssi=ssi, esi=esi,
                 csi=psi, opt=opt)
end
