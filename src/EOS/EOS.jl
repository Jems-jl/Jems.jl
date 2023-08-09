module EOS

using Jems.Chem, Jems.Constants

export AbstractEOS

abstract type AbstractEOS end

include("IdealEOS.jl")

end # module EOS
