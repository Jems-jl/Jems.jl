using LinearAlgebra

"""
    BicubicInterpolation <: AbstractBiInterpolation

Structure holding pre-calculated coefficients for bicubic (Catmull-Rom) interpolation.

# Fields
- `grid_x::Vector{Float64}`: Rectilinear grid axes for X.
- `grid_y::Vector{Float64}`: Rectilinear grid axes for Y.
- `coeff::Array{Float64, 4}`: 4D array `(num_vars, 16, Nx, Ny)` containing the 16 polynomial coefficients.
"""
struct BicubicInterpolation <: AbstractBiInterpolation
    grid_x :: Vector{Float64}
    grid_y :: Vector{Float64}
    coeff :: Array{Float64, 4} # (num_variables, 16, num_boxes_X, num_boxes_Y)

end


"""
    build_bicubic_interpolator(grid_x, grid_y, raw_data)

Calculates the 16 bicubic coefficients per box using a 4x4 stencil and Catmull-Rom logic.
Clamps indices at boundaries to prevent out-of-bounds errors. Uses `LinearAlgebra.mul!` 
for memory-efficient matrix multiplication. Returns an instance of `BicubicInterpolation`.
"""
function build_bicubic_interpolator(grid_x, grid_y, raw_data)
    num_vars, Nx_total, Ny_total = size(raw_data)
    Nx, Ny = length(grid_x) - 1, length(grid_y) - 1
    
    coeffs = Array{Float64, 4}(undef, num_vars, 16, Nx, Ny)
    P = Matrix{Float64}(undef, 4, 4)
    tmp = Matrix{Float64}(undef, 4, 4)
    C_mat = Matrix{Float64}(undef, 4, 4)
    
    M = 0.5 * [ 0.0  2.0  0.0  0.0;
               -1.0  0.0  1.0  0.0;
                2.0 -5.0  4.0 -1.0;
               -1.0  3.0 -3.0  1.0]
    @inbounds for var in 1:num_vars
        for j in 1:Ny 
            for i in 1:Nx 
                for n in 1:4 
                    for m in 1:4 
                        ix = clamp(i + m - 2, 1, Nx_total) # clamping any ouside index value to the boundary of the grid 
                        iy = clamp(j + n - 2, 1, Ny_total)
                        P[m, n] = raw_data[var, ix, iy]
                    end 
                end 
                mul!(tmp, M, P)
                mul!(C_mat, tmp, M')
                coeffs[var, :, i, j] .= vec(C_mat)
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
    C = view(interp.coeff, 1, :, i, j)
    # Horner's method for efficiency
    val1 = C[1] + u*(C[2] + u*(C[3] + u*C[4]))
    val2 = C[5] + u*(C[6] + u*(C[7] + u*C[8]))
    val3 = C[9] + u*(C[10] + u*(C[11] + u*C[12]))
    val4 = C[13] + u*(C[14] + u*(C[15] + u*C[16]))
    return val1 + v*(val2 + v*(val3 + v*val4))
end