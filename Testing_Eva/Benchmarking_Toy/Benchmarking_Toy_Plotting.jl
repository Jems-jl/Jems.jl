using BenchmarkTools
using Makie
using CairoMakie
using LinearAlgebra
using Statistics


##

min_v = 5
max_v = 50

##

ISO_files_A = []
ISO_files_B = []

file_indices = min_v:5:max_v

for i in file_indices
    push!(ISO_files_A, BenchmarkTools.load("BM1_Toy_$(i)ISO_A.json")[1])
    push!(ISO_files_B, BenchmarkTools.load("BM1_Toy_$(i)ISO_B.json")[1])
end

y_values_A = [median(ISO).time for ISO in ISO_files_A]
y_values_B = [median(ISO).time for ISO in ISO_files_B]

y_errors_A = [std(ISO).time ./ median(ISO).time for ISO in ISO_files_A]
y_errors_B = [std(ISO).time ./ median(ISO).time for ISO in ISO_files_B]

x_values = [Float64(i) for i in (min_v+4):5:(max_v+4)]

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


α, cov_matrix_A = weighted_least_squares_fit(log10.(x_values), log10.(y_values_A), y_errors_A)
β, cov_matrix_B = weighted_least_squares_fit(log10.(x_values), log10.(y_values_B), y_errors_B)

println("Intercept A: ", α[1])
println("Slope A: ", α[2])
println("Intercept B: ", β[1])
println("Slope B: ", β[2])


f = Figure();
ax = Axis(f[1, 1]; xlabel=L"N_{vars}", ylabel="Time")

α, cov_matrix_A = weighted_least_squares_fit(log10.(x_values), log10.(y_values_A), y_errors_A)
β, cov_matrix_B = weighted_least_squares_fit(log10.(x_values), log10.(y_values_B), y_errors_B)

fitted_line_A(x) = α[1] + α[2] * x
fitted_line_B(x) = β[1] + β[2] * x

# Plotting

errorbars!(log10.(x_values), log10.(y_values_A), y_errors_A, color=:blue, linewidth=3, whiskerwidth=10)
errorbars!(log10.(x_values), log10.(y_values_B), y_errors_B, color=:orange, linewidth=3, whiskerwidth=10)

scatter!(log10.(x_values), log10.(y_values_A), color=:blue, markersize = 15, marker = :diamond, label="1st benchmark")
scatter!(log10.(x_values), log10.(y_values_B), color=:orange, markersize = 15, marker = :diamond, label="2nd benchmark")

lines!(log10.(x_values), x -> fitted_line_A(x), color=:grey, linewidth=2, linestyle=:solid, label="Weighted Least Squares Fit")
lines!(log10.(x_values), x -> fitted_line_B(x), color=:grey, linewidth=2, linestyle=:solid)

axislegend(position=:lt)

# save("Benchmark_VaryingIsotopes_Toy_5to50ISO.png", f)
# save("Benchmark_VaryingIsotopes_Toy_20to50ISO.png", f)
# save("Benchmark_VaryingIsotopes_Toy_5to25ISO.png", f)
# save("Benchmark_VaryingIsotopes_Toy_5to20ISO.png", f)

f


##

