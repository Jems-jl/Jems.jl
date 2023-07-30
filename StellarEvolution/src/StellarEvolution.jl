module StellarEvolution

export StellarModel

using StellarConstants
using StellarChem
using StellarEOS
using StellarOpacity

using BlockBandedMatrices
using ForwardDiff
using Base.Threads

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
    mstar::Real # total model mass
    
    # Some basic info
    eos::StellarEOS.AbstractEOS
    opacity::StellarOpacity.AbstractOpacity
    isotope_data::Dict{Symbol, Isotope}

    # Jacobian matrix
    jac::BlockBandedMatrix
    function StellarModel(varnames::Vector{Symbol}, structure_equations::Vector{Function}, nvars::Int, nspecies::Int, nz, eos, opacity)
        ind_vars = ones(nvars*nz)
        eqs = ones(nvars*nz)
        m = ones(nz)

        l,u = 1,1          # block bandwidths
        N = M = nz        # number of row/column blocks
        cols = rows = [nvars for i in 1:N]  # block sizes

        jac = BlockBandedMatrix(Zeros(sum(rows),sum(cols)), rows,cols, (l,u))

        isotope_data = StellarChem.get_isotope_list();

        vari::Dict{Symbol,Int} = Dict()
        for i in eachindex(varnames)
            vari[varnames[i]] = i
        end

        new(ind_vars, varnames, eqs, nvars, nspecies, structure_equations, vari, nz, m, 0.0,eos,opacity,isotope_data,jac)
    end
end

include("Solver.jl")
include("Equations.jl")

end # module StellarModel
