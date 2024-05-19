✞ ☎
module Interp
# Code for interpolation for various orders
using LinearAlgebra
using Statistics
import Base.length
export CubicSpline, interp, slope, slope2, pchip, pchip2, pchip3
const eps = 1e-3 # rel error allowed on extrapolation
"""
CubicSpline(x,a,b,c,d)
concrete type for holding the data needed
to do a cubic spline interpolation
"""
abstract type AbstractSpline end
struct CubicSpline <: AbstractSpline
x::Union{Array{Float64,1},
StepRangeLen{Float64,
Base.TwicePrecision{Float64},
Base.TwicePrecision{Float64}}}
a::Array{Float64,1}
b::Array{Float64,1}
