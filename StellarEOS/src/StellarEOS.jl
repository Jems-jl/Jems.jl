module StellarEOS

using StellarChem, StellarConstants

abstract type AbstractEOS end

include("IdealEOS.jl")

end # module StellarEOS
