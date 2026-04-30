
"""
    BilinearInterpolation <: AbstractBiInterpolation

Structure holding pre-calculated coefficients for bilinear interpolation.

# Fields
- `grid_x::Vector{Float64}`: Rectilinear grid axes for X.
- `grid_y::Vector{Float64}`: Rectilinear grid axes for Y.
- `coeff::Array{Float64, 4}`: 4D array `(num_vars, 4, Nx, Ny)` containing linear coefficients `{C1, C2, C3, C4}`.
"""
struct BilinearInterpolation <: AbstractBiInterpolation
    grid_x :: Vector{Float64}
    grid_y :: Vector{Float64}
    coeff :: Array{Float64, 4} # (num_variables, 4, num_boxes_X, num_boxes_Y)
end 


"""
    build_bilinear_interpolator(grid_x, grid_y, raw_data)

Calculates the 4 bilinear coefficients for every rectangular box in the grid.
Assumes `raw_data` is dimensioned as `[var, x_index, y_index]`.Returns an 
instance of `BilinearInterpolation`.
"""
function build_bilinear_interpolator(grid_x, grid_y, raw_data) # raw_data[var, i, j].
    num_vars = size(raw_data, 1)
    Nx = length(grid_x) - 1 
    Ny = length(grid_y) - 1

    # Allocation for the coeff of interpolation 
    coeffs = Array{Float64, 4}(undef, num_vars, 4, Nx, Ny)

    for var in 1:num_vars
        for i in 1:Nx 
            for j in 1:Ny
                V00 = raw_data[var, i, j]
                V10 = raw_data[var, i+1, j]
                V01 = raw_data[var, i, j+1]
                V11 = raw_data[var, i+1, j+1]

                # f(x,y) = C1 + C2*x + C3*y + C4*x*y 
                C1 = V00 
                C2 = (V10 - V00)
                C3 = (V01 - V00) 
                C4 = (V11 - V01) - (V10 - V00)

                #storing the coefficient 
                coeffs[var,1, i, j] = C1
                coeffs[var,2, i, j] = C2
                coeffs[var,3, i, j] = C3
                coeffs[var,4, i, j] = C4
            end
        end 
    end 
    return BilinearInterpolation(grid_x, grid_y, coeffs)
end 

"""
    function evaluate_interp(interp, i, j, u, v)
Bilinear Evaluation function for the variable using the relative positioning in the box (u,v) and 
the interpolation coefficients for the same box 
"""

@inline function evaluate_interp(interp::BilinearInterpolation, i::Int, j::Int, u::T, v::T) where T
    C = view(interp.coeff, 1, :, i, j) #view just points instead of copying  
    # Polynomial: C1 + C2*u + C3*v + C4*u*v
    return C[1] + C[2]*u + C[3]*v + C[4]*u*v
end

