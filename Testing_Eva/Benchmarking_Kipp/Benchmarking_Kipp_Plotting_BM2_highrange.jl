using BenchmarkTools
using Makie
using CairoMakie
using LinearAlgebra
using Statistics


##

min_v = 15
max_v = 50
interval = 1

##

# Code to load in the first single version of the benchmarks

ISO_files_A = []
ISO_files_B = []

file_indices = min_v:interval:max_v

for i in file_indices
    push!(ISO_files_A, BenchmarkTools.load("BM1_Kipp_$(i)ISO_A.json")[1])
    push!(ISO_files_B, BenchmarkTools.load("BM1_Kipp_$(i)ISO_B.json")[1])
end

y_values_A = [minimum(ISO).time for ISO in ISO_files_A]
y_values_B = [minimum(ISO).time for ISO in ISO_files_B]

y_errors_A = [std(ISO).time ./ minimum(ISO).time for ISO in ISO_files_A]
y_errors_B = [std(ISO).time ./ minimum(ISO).time for ISO in ISO_files_B]

x_values = [Float64(i) for i in (min_v+4):interval:(max_v+4)]

##

# min_v = 15
# max_v = 50
# interval = 5

# ##

# # Code to load in the 2nd single version of the benchmarks V2 --> only every 5 isotopes

# ISO_files_A = []
# ISO_files_B = []

# file_indices = min_v:interval:max_v

# for i in file_indices
#     push!(ISO_files_A, BenchmarkTools.load("BM1_Kipp_$(i)ISO_A_V2.json")[1])
#     push!(ISO_files_B, BenchmarkTools.load("BM1_Kipp_$(i)ISO_B_V2.json")[1])
# end

# y_values_A_V2 = [minimum(ISO).time for ISO in ISO_files_A]
# y_values_B_V2 = [minimum(ISO).time for ISO in ISO_files_B]

# y_errors_A_V2 = [std(ISO).time ./ minimum(ISO).time for ISO in ISO_files_A]
# y_errors_B_V2 = [std(ISO).time ./ minimum(ISO).time for ISO in ISO_files_B]

# x_values_V2 = [Float64(i) for i in (min_v+4):interval:(max_v+4)]

##

# Redoing of the first datapoints
# The very first datapoints have quite a big error, so multiple benchmarks try to contribute
# to a smaller error. Only done for the first 4 datapoints


# min_v = 5
# max_v = 10
# interval = 1


# ISO_medians_A = []; ISO_medians_B = []
# ISO_errors_A = []; ISO_errors_B = []

# x_values_rep = []

# file_indices = min_v:interval:max_v

# for i in file_indices

#     ISO_file_A_ind = []
#     ISO_file_B_ind = []

#     push!(ISO_file_A_ind, BenchmarkTools.load("BM1_Kipp_$(i)ISO_A_V2_rep1.json")[1])
#     push!(ISO_file_A_ind, BenchmarkTools.load("BM1_Kipp_$(i)ISO_A_V2_rep2.json")[1])
#     push!(ISO_file_A_ind, BenchmarkTools.load("BM1_Kipp_$(i)ISO_A_V2_rep3.json")[1])

#     push!(ISO_file_B_ind, BenchmarkTools.load("BM1_Kipp_$(i)ISO_B_V2_rep1.json")[1])
#     push!(ISO_file_B_ind, BenchmarkTools.load("BM1_Kipp_$(i)ISO_B_V2_rep2.json")[1])
#     push!(ISO_file_B_ind, BenchmarkTools.load("BM1_Kipp_$(i)ISO_B_V2_rep3.json")[1])

#     y_medians_A_rep = [median(ISO).time for ISO in ISO_file_A_ind]
#     y_medians_B_rep = [median(ISO).time for ISO in ISO_file_B_ind]
    
#     y_gvalue_A_rep = [std(ISO).time^(-2) for ISO in ISO_file_A_ind]
#     y_gvalue_B_rep = [std(ISO).time^(-2) for ISO in ISO_file_B_ind]

#     sum_1_A = sum(x*y for (x, y) in zip((y_medians_A_rep), y_gvalue_A_rep))
#     sum_1_B = sum(x*y for (x, y) in zip((y_medians_B_rep), y_gvalue_B_rep))
    
#     sum_2_A = sum(y_gvalue_A_rep)
#     sum_2_B = sum(y_gvalue_B_rep)


#     # println("sum_1 is: ", sum_1)
#     # println("sum_2 is: ", sum_2)

#     m_A = sum_1_A / sum_2_A
#     m_B = sum_1_B / sum_2_B
#     ϵ_A = (sum_2_A)^(-1/2)
#     ϵ_B = (sum_2_B)^(-1/2)

#     x_value = i + 4

#     # println("m is: ", m)
#     # println("ϵ is: ", ϵ)

#     push!(ISO_medians_A, m_A)
#     push!(ISO_medians_B, m_B)

#     push!(ISO_errors_A, ϵ_A)
#     push!(ISO_errors_B, ϵ_B)

#     push!(x_values_rep, x_value)

# end

# y_values_rep_A = ISO_medians_A
# y_errors_rep_A = ISO_errors_A ./ ISO_medians_A

# y_values_rep_B = ISO_medians_B
# y_errors_rep_B = ISO_errors_B ./ ISO_medians_B

# # println("y values are: ", y_values_rep_A)
# # println("y errors are: ", y_errors_rep_A)

##

# Updating the first few datapoints

# y_values_A[1:6] .= y_values_rep_A
# y_errors_A[1:6] .= y_errors_rep_A

# y_values_B[1:6] .= y_values_rep_B
# y_errors_B[1:6] .= y_errors_rep_B


##


# Set up Makie

using CairoMakie, LaTeXStrings, MathTeXEngine
basic_theme = Theme(fonts=(regular=texfont(:text), bold=texfont(:bold),
                           italic=texfont(:italic), bold_italic=texfont(:bolditalic)),
                    fontsize=30, size=(1000, 750), linewidth=7,
                    Axis=(xlabelsize=40, ylabelsize=40, titlesize=40, xgridvisible=true, ygridvisible=true,
                          spinewidth=2.5, xminorticksvisible=true, yminorticksvisible=true, xtickalign=1, ytickalign=1,
                          xminortickalign=1, yminortickalign=1, xticksize=14, xtickwidth=2.5, yticksize=14,
                          ytickwidth=2.5, xminorticksize=7, xminortickwidth=2.5, yminorticksize=7, yminortickwidth=2.5,
                          xticklabelsize=35, yticklabelsize=35, xticksmirrored=true, yticksmirrored=true),
                    Legend=(patchsize=(70, 10), framevisible=false, patchlabelgap=20, rowgap=10))
set_theme!(basic_theme)

##


##

function weighted_least_squares_fit(xdata, ydata, yerrors)
    @assert length(xdata) == length(ydata) == length(yerrors)
    n = length(xdata)
    X = hcat(ones(n), xdata) 
    W = Diagonal(1 ./ (yerrors .^ 2))
    β = W * X \ (W * ydata)
    cov_matrix = inv(X' * W * X)
    return β, cov_matrix
end

# V1

α, cov_matrix_A = weighted_least_squares_fit(log10.(x_values), log10.(y_values_A), y_errors_A)
β, cov_matrix_B = weighted_least_squares_fit(log10.(x_values), log10.(y_values_B), y_errors_B)

println("Intercept A: ", α[1])
println("Slope A: ", α[2])

println("Intercept B: ", β[1])
println("Slope B: ", β[2])


f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(\mathrm{N_{var}})", ylabel=L"\log_{10}(\mathrm{t/[ns]})")

α, cov_matrix_A = weighted_least_squares_fit(log10.(x_values), log10.(y_values_A), y_errors_A)
β, cov_matrix_B = weighted_least_squares_fit(log10.(x_values), log10.(y_values_B), y_errors_B)

fitted_line_A(x) = α[1] + α[2] * x
fitted_line_B(x) = β[1] + β[2] * x

# V2

# α_V2, cov_matrix_A_V2 = weighted_least_squares_fit(log10.(x_values_V2), log10.(y_values_A_V2), y_errors_A_V2)
# β_V2, cov_matrix_B_V2 = weighted_least_squares_fit(log10.(x_values_V2), log10.(y_values_B_V2), y_errors_B_V2)

# println("Intercept A V2: ", α_V2[1])
# println("Slope A V2: ", α_V2[2])

# println("Intercept B V2: ", β_V2[1])
# println("Slope B V2: ", β_V2[2])


# f = Figure();
# ax = Axis(f[1, 1]; xlabel=L"N_{vars}", ylabel="Time")

# α_V2, cov_matrix_A_V2 = weighted_least_squares_fit(log10.(x_values_V2), log10.(y_values_A_V2), y_errors_A_V2)
# β_V2, cov_matrix_B_V2 = weighted_least_squares_fit(log10.(x_values_V2), log10.(y_values_B_V2), y_errors_B_V2)

# fitted_line_A_V2(x) = α_V2[1] + α_V2[2] * x
# fitted_line_B_V2(x) = β_V2[1] + β_V2[2] * x

# Plotting

# V1

# errorbars!(log10.(x_values), log10.(y_values_A), y_errors_A, color=:blue, linewidth=2, whiskerwidth=10)
errorbars!(log10.(x_values), log10.(y_values_B), y_errors_B, color=:darkorange, linewidth=2, whiskerwidth=10)

# scatter!(log10.(x_values), log10.(y_values_A), color=:blue, markersize = 15, marker = :diamond, label="Benchmark Jacobian evaluation")
scatter!(log10.(x_values), log10.(y_values_B), color=:darkorange, markersize = 15, marker = :diamond, label="Benchmark Jacobian + solver")

# lines!(log10.(x_values), x -> fitted_line_A(x), color=:grey, linewidth=2, linestyle=:solid, label="Fit")
lines!(log10.(x_values), x -> fitted_line_B(x), color=:grey, linewidth=2, linestyle=:solid, label="Fit")

# V2

# errorbars!(log10.(x_values_V2), log10.(y_values_A_V2), y_errors_A_V2, color=:red, linewidth=3, whiskerwidth=10)
# errorbars!(log10.(x_values_V2), log10.(y_values_B_V2), y_errors_B_V2, color=:green, linewidth=3, whiskerwidth=10)

# scatter!(log10.(x_values_V2), log10.(y_values_A_V2), color=:red, markersize = 15, marker = :diamond, label="1st benchmark V2")
# scatter!(log10.(x_values_V2), log10.(y_values_B_V2), color=:green, markersize = 15, marker = :diamond, label="2nd benchmark V2")

# lines!(log10.(x_values_V2), x -> fitted_line_A_V2(x), color=:black, linewidth=2, linestyle=:solid, label="Weighted Least Squares Fit")
# lines!(log10.(x_values_V2), x -> fitted_line_B_V2(x), color=:black, linewidth=2, linestyle=:solid)

# Repeated datapoints

# errorbars!(log10.(x_values_rep), log10.(y_values_rep_A), y_errors_rep_A, color=:red, linewidth=3, whiskerwidth=10)
# errorbars!(log10.(x_values_rep), log10.(y_values_rep_B), y_errors_rep_B, color=:green, linewidth=3, whiskerwidth=10)

# scatter!(log10.(x_values_rep), log10.(y_values_rep_A), color=:red, markersize = 15, marker = :diamond, label="1st benchmark")
# scatter!(log10.(x_values_rep), log10.(y_values_rep_B), color=:green, markersize = 15, marker = :diamond, label="2nd benchmark")


axislegend(position=:lt)

f


##

save("BM2_VaryingIsotopes_Kipp_15_to_50.png", f)

# Intercept A: 4.395836352203541
# Slope A: 2.3140012129831296
# Intercept B: 4.655071721517739
# Slope B: 2.2725071973635687