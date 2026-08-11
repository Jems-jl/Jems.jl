using LinearAlgebra

"""
    BicubicInterpolation <: Abstract2DInterpolation

Structure holding pre-calculated coefficients for bicubic (Catmull-Rom) interpolation.

# Fields
- `grid_x::Vector{Float64}`: Rectilinear grid axes for X.
- `grid_y::Vector{Float64}`: Rectilinear grid axes for Y.
- `coeff::Array{Float64, 4}`: 4D array `(num_vars, 16, Nx, Ny)` containing the 16 polynomial coefficients.
"""
struct BicubicInterpolation <: Abstract2DInterpolation
    grid_x :: Vector{Float64}
    grid_y :: Vector{Float64}
    coeff :: Array{Float64, 3} # (num_variables, 16, num_boxes_X, num_boxes_Y)

end


"""
    build_bicubic_interpolator(grid_x, grid_y, raw_data)

Calculates the 16 bicubic coefficients per box using a 4x4 stencil and Catmull-Rom logic.
Clamps indices at boundaries to prevent out-of-bounds errors. Uses `LinearAlgebra.mul!` 
for memory-efficient matrix multiplication. Returns an instance of `BicubicInterpolation`.
"""
function build_bicubic_interpolator(grid_x, grid_y, raw_data::AbstractArray{<:Real, 3})
    num_vars, Nx_total, Ny_total = size(raw_data)
    Nx, Ny = length(grid_x) - 1, length(grid_y) - 1
    
    coeffs = Array{Float64, 3}(undef, 16, Nx, Ny)
    P = Matrix{Float64}(undef, 4, 4)
    tmp = Matrix{Float64}(undef, 4, 4)
    C_mat = Matrix{Float64}(undef, 4, 4)
    
    M = 0.5 * [ 0.0  2.0  0.0  0.0;
               -1.0  0.0  1.0  0.0;
                2.0 -5.0  4.0 -1.0;
               -1.0  3.0 -3.0  1.0]
               
    Mt = copy(M') # Explicitly precompute the adjoint transpose outside the loop
    
    @inbounds for j in 1:Ny 
        for i in 1:Nx 
            for n in 1:4 
                for m in 1:4 
                    ix = clamp(i + m - 2, 1, Nx_total) 
                    iy = clamp(j + n - 2, 1, Ny_total)
                    P[m, n] = raw_data[1, ix, iy] # Var is always 1 because we pass a slice
                end 
            end 
            mul!(tmp, M, P)
            mul!(C_mat, tmp, Mt)
            
            # Zero-allocation manual flattening (replaces .= vec(C_mat) slicing)
            idx = 1
            for cn in 1:4
                for cm in 1:4
                    coeffs[idx, i, j] = C_mat[cm, cn]
                    idx += 1
                end
            end
        end 
    end 
    return BicubicInterpolation(grid_x, grid_y, coeffs)
end
"""
    function evaluate_interp(interp, i, j, u, v)

Evaluates the bicubic polynomial for a specific variable `var` at fractional coordinates `(u, v)` 
inside the grid box `(i, j)`. Uses Horner's method for better execution speed.
"""

@inline function evaluate_interp(interp::BicubicInterpolation, i::Int, j::Int, u::T, v::T) where T
    # Explicit unpacking is infinitely faster than a runtime view()
    @inbounds begin
        c1  = interp.coeff[1, i, j]
        c2  = interp.coeff[2, i, j]
        c3  = interp.coeff[3, i, j]
        c4  = interp.coeff[4, i, j]
        c5  = interp.coeff[5, i, j]
        c6  = interp.coeff[6, i, j]
        c7  = interp.coeff[7, i, j]
        c8  = interp.coeff[8, i, j]
        c9  = interp.coeff[9, i, j]
        c10 = interp.coeff[10, i, j]
        c11 = interp.coeff[11, i, j]
        c12 = interp.coeff[12, i, j]
        c13 = interp.coeff[13, i, j]
        c14 = interp.coeff[14, i, j]
        c15 = interp.coeff[15, i, j]
        c16 = interp.coeff[16, i, j]
    end

    # Horner's method
    val1 = c1  + u*(c2  + u*(c3  + u*c4 ))
    val2 = c5  + u*(c6  + u*(c7  + u*c8 ))
    val3 = c9  + u*(c10 + u*(c11 + u*c12))
    val4 = c13 + u*(c14 + u*(c15 + u*c16))
    return val1 + v*(val2 + v*(val3 + v*val4))
end