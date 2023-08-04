using BlockBandedMatrices
using SparseArrays
using LinearSolve

mutable struct StellarStepInfo
    dt::Real
    ind_vars_old::Vector{<:Real}
    lnT_old::Vector{<:Real}
    lnL_old::Vector{<:Real}
    lnP_old::Vector{<:Real}
    lnρ_old::Vector{<:Real}
    lnr_old::Vector{<:Real}
end

"""
    StellarModel

Structured type that defines the state of the stellar model 
"""
mutable struct StellarModel
    # properties that define the model
    ind_vars::Vector{<:Real}
    varnames::Vector{Symbol}
    eqs::Vector{<:Real}
    nvars::Int # this is the sum of hydro vars and species
    nspecies::Int # Just the number of species in the network
    structure_equations::Vector{Function} # Be careful, this probably kills type inference
    vari::Dict{Symbol,Int} #links variable names to ind_vars vector

    # grid properties
    nz::Int # number of zones in the model
    m::Vector{<:Real} # mass coordinate of each cell
    dm::Vector{<:Real} # mass contained in each cell
    mstar::Real # total model mass
    
    # Some basic info
    eos::StellarEOS.AbstractEOS
    opacity::StellarOpacity.AbstractOpacity
    isotope_data::Dict{Symbol, Isotope}

    # Jacobian matrix
    jac::SparseMatrixCSC{Float64, Int64}
    linear_solver  #solver that is produced by LinearSolve

    # Information computed at the start of the Step
    ssi::StellarStepInfo

    # Space for used defined options, defaults are in Options.jl
    opt::Options
    function StellarModel(varnames::Vector{Symbol}, structure_equations::Vector{Function}, nvars::Int, nspecies::Int, nz, eos, opacity)
        ind_vars = ones(nvars*nz)
        eqs = ones(nvars*nz)
        m = ones(nz)

        l,u = 1,1          # block bandwidths
        N = M = nz        # number of row/column blocks
        cols = rows = [nvars for i in 1:N]  # block sizes

        jac_BBM = BlockBandedMatrix(Ones(sum(rows),sum(cols)), rows,cols, (l,u))
        jac = sparse(jac_BBM)
        #create solver
        problem = LinearProblem(jac, eqs)
        linear_solver = init(problem)

        isotope_data = StellarChem.get_isotope_list();

        vari::Dict{Symbol,Int} = Dict()
        for i in eachindex(varnames)
            vari[varnames[i]] = i
        end

        dm = zeros(nz)

        ssi = StellarStepInfo(0.0,[0.0],[0.0],[0.0],[0.0],[0.0],[0.0])

        opt = Options()

        new(ind_vars, varnames, eqs, nvars, nspecies, structure_equations, vari, nz, m, dm, 0.0,
            eos,opacity,isotope_data,jac,linear_solver,ssi, opt)
    end
end
