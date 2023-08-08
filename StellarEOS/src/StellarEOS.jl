module StellarEOS

using StellarChem, StellarConstants

export AbstractEOS

abstract type AbstractEOS end

include("IdealEOS.jl")

end # module StellarEOS
