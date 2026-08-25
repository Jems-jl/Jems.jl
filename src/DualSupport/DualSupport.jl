module DualSupport

using ForwardDiff
using StaticArrays

export LocalDualData, update_local_dual_data_value!, update_local_dual_data!,
        get_local_dual, get_m1_dual, get_00_dual, get_p1_dual, get_value, 
        get_slope, construct_derivative, interpolate_to_value

# Inspired by DiffCache from PreallocationTools (https://github.com/SciML/PreallocationTools.jl)
"""
    struct StarDiffCache{SIZE, TNUMBER, DUAL_TAG}

Definition of StarDiffCache, a cache that makes room to store partial derivatives.
Parametric in types `SIZE`, the size of the array, and `TNUMBER`, the type of the number used for calculations. 
DUAL_TAG is the tag used to identify the Dual number order, required to use nested duals.
"""
struct StarDiffCache{SIZE,TNUMBER,DUAL_TAG}
    dual_data::MVector{SIZE,TNUMBER}
    dual_tag::Type{DUAL_TAG}
end


"""
    function StarDiffCache(nvars::Int, ::Type{TNUMBER}, ::TYPE{DUAL_TAG}) where {TNUMBER}

Instantiates a StarDiffCache object of size `nvars+1`, and fills it with zeros.
"""
function StarDiffCache(nvars::Int, ::Type{TNUMBER}, ::Type{DUAL_TAG}) where {TNUMBER,DUAL_TAG}
    StarDiffCache{nvars + 1,TNUMBER,DUAL_TAG}(zeros(TNUMBER, nvars + 1), DUAL_TAG)
end

## This uses reinterpret
#function get_dual(sdc::StarDiffCache{SIZE, TNUMBER}) where{SIZE,TNUMBER}
#    reinterpret(ForwardDiff.Dual{Nothing, TNUMBER, SIZE}, sdc.dual_data)[1]
#end

# kudos to user Mason Protter from discourse.julia.com
# beware of caveats
# https://discourse.julialang.org/t/reinterpret-vector-into-single-struct/107709
function get_dual(sdc::StarDiffCache{SIZE,TNUMBER,DUAL_TAG}) where {SIZE,TNUMBER,DUAL_TAG}
    p::Ptr{ForwardDiff.Dual{DUAL_TAG,TNUMBER,SIZE-1}} = pointer(sdc.dual_data)
    unsafe_load(p)         # Load the first element from that pointer
end

function get_mixed_dual(sdc::StarDiffCache{SIZE,TNUMBER,DUAL_TAG}) where {SIZE,TNUMBER,DUAL_TAG}
    p::Ptr{ForwardDiff.Dual{DUAL_TAG,TNUMBER,(SIZE-1)*2÷3}} = pointer(sdc.dual_data)
    unsafe_load(p)         # Load the first element from that pointer
end

"""
Interpolate the linear slope of y to x at x[i]
"""
function get_slope(xs, ys, i)
    if i == length(xs)
        throw(ArgumentError("i=$i is the last index of xs, cannot compute slope"))
    end
    return (ys[i+1] - ys[i]) / (xs[i+1] - xs[i])
end

function interpolate_to_value(x1, x2, y1, y2, x)
    if x < min(ForwardDiff.value(x1), ForwardDiff.value(x2)) || x > max(ForwardDiff.value(x1), ForwardDiff.value(x2))  # coerce value here because we don't want partials to affect the bounds check
        throw(ArgumentError("x=$x is outside the range [$x1, $x2]"))
    end
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
end


"""
Compute the value 
    dy/dx|_{z} = ∂y/∂x - Δy/Δz * ∂z/∂x
per the implicit function theorem.
"""
function construct_derivative(partial_y, partial_z, slope)
    return partial_y - slope * partial_z
end


include("LocalDualData.jl")
include("MixedDualData.jl")

end
