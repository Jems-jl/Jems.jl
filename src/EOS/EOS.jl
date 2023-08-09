module EOS

using ..Chem, ..Constants

export AbstractEOS

abstract type AbstractEOS end

include("IdealEOS.jl")

end # module EOS
