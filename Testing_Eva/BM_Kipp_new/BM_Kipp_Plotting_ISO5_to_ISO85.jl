using LinearAlgebra
using BenchmarkTools
using Statistics

n_list = [];  
m_list_A = []; m_list_B = [];  
s_list_A = []; s_list_B = [];
q10_list_A = []; q10_list_B = [];
q25_list_A = []; q25_list_B = [];
q75_list_A = []; q75_list_B = [];
q90_list_A = []; q90_list_B = [];

interval_1 = 5:1:20
interval_2 = 25:5:120

##

# reload the values from the benchmarks 5:1:20 and 22:2:100

for n in interval_1

    bm_A = BenchmarkTools.load("BM1_Kipp_$(n)ISO_A_V2.json")[1]
    bm_B = BenchmarkTools.load("BM1_Kipp_$(n)ISO_B_V2.json")[1]

    data_A = bm_A.times
    data_B = bm_B.times

    m_A = median(data_A); m_B = median(data_B)
    # println("median is: ", m)
    perc_10_A = quantile(data_A, 0.10); perc_10_B = quantile(data_B, 0.10)
    perc_25_A = quantile(data_A, 0.25); perc_25_B = quantile(data_B, 0.25)
    perc_75_A = quantile(data_A, 0.75); perc_75_B = quantile(data_B, 0.75)
    perc_90_A = quantile(data_A, 0.90); perc_90_B = quantile(data_B, 0.90)

    push!(n_list, n+4)

    push!(m_list_A, m_A); push!(m_list_B, m_B)

    push!(q10_list_A, perc_10_A); push!(q10_list_B, perc_10_B)
    push!(q25_list_A, perc_25_A); push!(q25_list_B, perc_25_B)
    push!(q75_list_A, perc_75_A); push!(q75_list_B, perc_75_B)
    push!(q90_list_A, perc_90_A); push!(q90_list_B, perc_90_B)


end


for n in interval_2

    bm_A = BenchmarkTools.load("BM1_Kipp_$(n)ISO_A_V2.json")[1]
    bm_B = BenchmarkTools.load("BM1_Kipp_$(n)ISO_B_V2.json")[1]

    data_A = bm_A.times
    data_B = bm_B.times

    m_A = median(data_A); m_B = median(data_B)
    # println("median is: ", m)
    perc_10_A = quantile(data_A, 0.10); perc_10_B = quantile(data_B, 0.10)
    perc_25_A = quantile(data_A, 0.25); perc_25_B = quantile(data_B, 0.25)
    perc_75_A = quantile(data_A, 0.75); perc_75_B = quantile(data_B, 0.75)
    perc_90_A = quantile(data_A, 0.90); perc_90_B = quantile(data_B, 0.90)

    push!(n_list, n+4)

    push!(m_list_A, m_A); push!(m_list_B, m_B)

    push!(q10_list_A, perc_10_A); push!(q10_list_B, perc_10_B)
    push!(q25_list_A, perc_25_A); push!(q25_list_B, perc_25_B)
    push!(q75_list_A, perc_75_A); push!(q75_list_B, perc_75_B)
    push!(q90_list_A, perc_90_A); push!(q90_list_B, perc_90_B)


end


##

y_err_top_A = m_list_A .- q25_list_A 
y_err_bot_A = q75_list_A .- m_list_A

y_err_top_B = m_list_B .- q25_list_B 
y_err_bot_B = q75_list_B .- m_list_B

##

x_values = n_list
y_values_A = m_list_A
y_values_B = m_list_B
# y_errors = s_list ./ m_list
y_errors_top_A = y_err_top_A ./ m_list_A
y_errors_bot_A = y_err_bot_A ./ m_list_A
y_errors_top_B = y_err_top_B ./ m_list_B
y_errors_bot_B = y_err_bot_B ./ m_list_B

##

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

# Fits maken --> niet relevant voor de grote code

# function weighted_least_squares_fit(xdata, ydata, yerrors_top, yerrors_bottom)
#     @assert length(xdata) == length(ydata) == length(yerrors_top) == length(yerrors_bottom)
#     n = length(xdata)
#     X = hcat(ones(n), xdata) 
#     W = Diagonal(1 ./ ((yerrors_top .^ 2) + (yerrors_bottom .^ 2)))  # Combined weights
#     ydata_weighted = (yerrors_top .* ydata .+ yerrors_bottom .* ydata) ./ (yerrors_top .+ yerrors_bottom)  # Weighted ydata
#     β = W * X \ (W * ydata_weighted)
#     cov_matrix = inv(X' * W * X)
#     return β, cov_matrix
# end

# function fit_errors(cov_matrix)
#     errors = sqrt.(diag(cov_matrix))
#     return errors
# end

# α, cov_matrix = weighted_least_squares_fit(log10.(x_values), log10.(y_values), y_errors_top, y_errors_bot)

# println("Intercept A: ", α[1])
# println("Slope A: ", α[2])
# println("Error on Intercept A: ", fit_errors(cov_matrix)[1])
# println("Error on Slope A: ", fit_errors(cov_matrix)[2])

# fitted_line_A(x) = α[1] + α[2] * x

##

# Fits uit de andere files hier ingeladen en geplot

fitted_line_1_A(x) = 5.510 + 1.453 * x
fitted_line_1_B(x) = 5.339 + 1.737 * x


# Intercept A: 5.510518627479121
# Slope A: 1.453384603037347
# Error on Intercept A: 0.20620990484447807
# Error on Slope A: 0.1763784448452581
# Intercept B: 5.339619971479433
# Slope B: 1.7374830280945186
# Error on Intercept B: 0.22289229775096414
# Error on Slope B: 0.19420386966801031


fitted_line_2_A(x) = 4.194 + 2.406 * x
fitted_line_2_B(x) = 4.298 + 2.511 * x

# Intercept A: 4.194749393645936
# Slope A: 2.406600625009459
# Error on Intercept A: 0.08046371818924657
# Error on Slope A: 0.047524783945580296
# Intercept B: 4.29811958609154
# Slope B: 2.511842075504802
# Error on Intercept B: 0.12488946540351589
# Error on Slope B: 0.07491207924636788

##

# Aangepaste x-ranges voor de verschillende plots

x_values_1 = [x for x in x_values if x < 30]
x_values_2 = [x for x in x_values if x > 20]
# x_values_3 = [x for x in x_values if x > 55]

##

f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(\mathrm{n})", ylabel=L"\log_{10}(\mathrm{t/[ns]})")

lines!(log10.(x_values_1), x -> fitted_line_1_A(x), color=:grey, linewidth=2, linestyle=:solid, label="Fit")
lines!(log10.(x_values_2), x -> fitted_line_2_A(x), color=:grey, linewidth=2, linestyle=:solid, label="Fit")

lines!(log10.(x_values_1), x -> fitted_line_1_B(x), color=:grey, linewidth=2, linestyle=:solid, label="Fit")
lines!(log10.(x_values_2), x -> fitted_line_2_B(x), color=:grey, linewidth=2, linestyle=:solid, label="Fit")

errorbars!(log10.(x_values), log10.(y_values_A), y_errors_top_A, y_errors_bot_A, color=:blue, linewidth=2, whiskerwidth=10)
scatter!(log10.(x_values), log10.(y_values_A), color=:blue, markersize = 15, marker = :diamond, label="Benchmark Jacobian + solver")

errorbars!(log10.(x_values), log10.(y_values_B), y_errors_top_B, y_errors_bot_B, color=:orange, linewidth=2, whiskerwidth=10)
scatter!(log10.(x_values), log10.(y_values_B), color=:orange, markersize = 15, marker = :diamond, label="Benchmark Jacobian + solver")


# ylims!(ax, 1.0, 5.0)

f

##

# save("BM_MM_n5_to_n100_fit.png", f)

##

