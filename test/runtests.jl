using Jems
using Jems.Interpolations
using Test

@testset "Jems.jl" begin
    @testset "Bilinear Interpolation" begin
        grid_x = [0.0, 1.0, 2.0]
        grid_y = [0.0, 1.0]
        data = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], (1, 3, 2))
        
        interp = build_bilinear_interpolator(grid_x, grid_y, data)
        
        i, j, u, v = get_data_position(grid_x, grid_y, 0.5, 0.3)
        result = evaluate_interp(interp, i, j, u, v)
        
        @test result ≈ 2.4 rtol=1e-2
    end
    
    @testset "Bicubic Interpolation" begin

        grid_x = [0.0, 1.0, 2.0, 3.0, 4.0]
        grid_y = [0.0, 1.0, 2.0, 3.0]
        
       
        data = zeros(1, 5, 4)
        for i in 1:5
            for j in 1:4
                data[1, i, j] = (grid_x[i]^2 + grid_y[j]^2)
            end
        end
        
        interp = build_bicubic_interpolator(grid_x, grid_y, data)
        
        i, j, u, v = get_data_position(grid_x, grid_y, 1.5, 1.5)
        result = evaluate_interp(interp, i, j, u, v)
        
        @test result ≈ 4.5 rtol=1e-1  
    end
end
