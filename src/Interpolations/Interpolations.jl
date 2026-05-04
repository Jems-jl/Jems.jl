module Interpolations

using ForwardDiff
using LinearAlgebra

export Abstract2DInterpolation, get_data_position, 
       BilinearInterpolation, build_bilinear_interpolator,
       BicubicInterpolation, build_bicubic_interpolator, evaluate_interp

"""
    abstract type Abstract2DInterpolation

Abstract supertype from which all Bi-Interpolation definitions must derive

`struct BilinearInterpolation <: Abstract2DInterpolation`
"""
abstract type Abstract2DInterpolation end

"""
    get_data_position(grid_x, grid_y, x, y)

Locates the box indices (i, j) and calculates fractional offsets (u, v) in [0, 1],
given the x,y values in the (grid_x, grid_y)
"""

function get_data_position(grid_x, grid_y, x::T, y::T) where T
    i = clamp(searchsortedlast(grid_x, ForwardDiff.value(x)), 1, length(grid_x) - 1)
    j = clamp(searchsortedlast(grid_y, ForwardDiff.value(y)), 1, length(grid_y) - 1)

    # fractional position in the data stencil 
    u = (x - grid_x[i])/ (grid_x[i+1]-grid_x[i])
    v = (y - grid_y[j]) / (grid_y[j+1] - grid_y[j])

    return i, j, u, v
end 

include("BilinearInterpolation.jl")
include("BicubicInterpolation.jl")

end 