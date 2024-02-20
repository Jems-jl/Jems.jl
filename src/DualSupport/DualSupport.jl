module DualSupport

using ForwardDiff
using StaticArrays

# Inspired by DiffCache from PreallocationTools (https://github.com/SciML/PreallocationTools.jl)
struct StarDiffCache{SIZE, TNUMBER}
    dual_data::MVector{SIZE,TNUMBER}
end

function StarDiffCache(nvars::Int, ::Type{TNUMBER}) where {TNUMBER}
    StarDiffCache{nvars+1, TNUMBER}(zeros(TNUMBER, nvars+1))
end

## This uses reinterpret
#function get_dual(sdc::StarDiffCache{SIZE, TNUMBER}) where{SIZE,TNUMBER}
#    reinterpret(ForwardDiff.Dual{Nothing, TNUMBER, SIZE}, sdc.dual_data)[1]
#end

# kudos to user Mason Protter from discourse.julia.com
# beware of caveats
# https://discourse.julialang.org/t/reinterpret-vector-into-single-struct/107709
function get_dual(sdc::StarDiffCache{SIZE, TNUMBER}) where {SIZE,TNUMBER}
    p::Ptr{ForwardDiff.Dual{Nothing, TNUMBER, SIZE-1}} = pointer(sdc.dual_data)
    unsafe_load(p)         # Load the first element from that pointer
end

function get_face_dual(sdc::StarDiffCache{SIZE, TNUMBER}) where {SIZE,TNUMBER}
    p::Ptr{ForwardDiff.Dual{Nothing, TNUMBER, (SIZE-1)*2÷3}} = pointer(sdc.dual_data)
    unsafe_load(p)         # Load the first element from that pointer
end

include("CellDualData.jl")
include("FaceDualData.jl")

end