using LinearAlgebra
using BenchmarkTools
using Statistics

n_list = [];  # x-values of n
m_list = [];  # medians of the Benchmarks
s_list = [];
q10_list = [];
q25_list = [];
q75_list = [];
q90_list = [];

interval = 62:2:100

##

# code for a single BM per datapoint

# for i in interval

#     n = 10^(i)

#     A = rand(n,n); B = rand(n,n); C = rand(n,n)
#     bm = @benchmark mul!($A,$B,$C)
#     BenchmarkTools.save("BM_MM_n$(n).json", bm)

#     data = bm.times
#     s = std(bm).time

#     m = median(data)
#     println("median is: ", m)
#     perc_10 = quantile(data, 0.10)
#     perc_25 = quantile(data, 0.25)
#     perc_75 = quantile(data, 0.75)
#     perc_90 = quantile(data, 0.90)


#     push!(n_list, n)
#     push!(m_list, m)
#     push!(s_list, s)
#     push!(q10_list, perc_10)
#     push!(q25_list, perc_25)
#     push!(q75_list, perc_75)
#     push!(q90_list, perc_90)

# end

##

# code for multiple takes of a BM

for n in interval

    A = rand(n,n); B = rand(n,n); C = rand(n,n)
    bm1 = @benchmark mul!($A,$B,$C)
    BenchmarkTools.save("BM_MM_n$(n)_v1.json", bm1)

    A = rand(n,n); B = rand(n,n); C = rand(n,n)
    bm2 = @benchmark mul!($A,$B,$C)
    BenchmarkTools.save("BM_MM_n$(n)_v2.json", bm1)

    A = rand(n,n); B = rand(n,n); C = rand(n,n)
    bm3 = @benchmark mul!($A,$B,$C)
    BenchmarkTools.save("BM_MM_n$(n)_v3.json", bm1)
    
    data = vcat(bm1.times, bm2.times, bm3.times)

    m = median(data)
    println("median is: ", m)
    perc_10 = quantile(data, 0.10)
    perc_25 = quantile(data, 0.25)
    perc_75 = quantile(data, 0.75)
    perc_90 = quantile(data, 0.90)

    push!(n_list, n)
    push!(m_list, m)
    # push!(s_list, s)
    push!(q10_list, perc_10)
    push!(q25_list, perc_25)
    push!(q75_list, perc_75)
    push!(q90_list, perc_90)

end

##

# reload the values from the benchmarks 62:2:100

# file_indices = 62:2:100

# for n in file_indices

#     bm1 = BenchmarkTools.load("BM_MM_n$(n)_v1.json")[1]
#     bm2 = BenchmarkTools.load("BM_MM_n$(n)_v2.json")[1]
#     bm3 = BenchmarkTools.load("BM_MM_n$(n)_v3.json")[1]

#     data = vcat(bm1.times, bm2.times, bm3.times)

#     m = median(data)
#     # println("median is: ", m)
#     perc_10 = quantile(data, 0.10)
#     perc_25 = quantile(data, 0.25)
#     perc_75 = quantile(data, 0.75)
#     perc_90 = quantile(data, 0.90)

#     push!(n_list, n)
#     push!(m_list, m)
#     # push!(s_list, s)

#     push!(q10_list, perc_10)
#     push!(q25_list, perc_25)
#     push!(q75_list, perc_75)
#     push!(q90_list, perc_90)


# end



##

y_err_top = m_list .- q25_list 
y_err_bot = q75_list .- m_list

##

x_values = n_list
y_values = m_list
# y_errors = s_list ./ m_list
y_errors_top = y_err_top ./ m_list
y_errors_bot = y_err_bot ./ m_list

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

function weighted_least_squares_fit(xdata, ydata, yerrors_top, yerrors_bottom)
    @assert length(xdata) == length(ydata) == length(yerrors_top) == length(yerrors_bottom)
    n = length(xdata)
    X = hcat(ones(n), xdata) 
    W = Diagonal(1 ./ ((yerrors_top .^ 2) + (yerrors_bottom .^ 2)))  # Combined weights
    ydata_weighted = (yerrors_top .* ydata .+ yerrors_bottom .* ydata) ./ (yerrors_top .+ yerrors_bottom)  # Weighted ydata
    β = W * X \ (W * ydata_weighted)
    cov_matrix = inv(X' * W * X)
    return β, cov_matrix
end

function fit_errors(cov_matrix)
    errors = sqrt.(diag(cov_matrix))
    return errors
end

α, cov_matrix = weighted_least_squares_fit(log10.(x_values), log10.(y_values), y_errors_top, y_errors_bot)

println("Intercept A: ", α[1])
println("Slope A: ", α[2])
println("Error on Intercept A: ", fit_errors(cov_matrix)[1])
println("Error on Slope A: ", fit_errors(cov_matrix)[2])

fitted_line_A(x) = α[1] + α[2] * x

##

f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(\mathrm{n})", ylabel=L"\log_{10}(\mathrm{t/[ns]})")

lines!(log10.(x_values), x -> fitted_line_A(x), color=:grey, linewidth=2, linestyle=:solid, label="Fit")
errorbars!(log10.(x_values), log10.(y_values), y_errors_top, y_errors_bot, color=:blue, linewidth=2, whiskerwidth=10)
scatter!(log10.(x_values), log10.(y_values), color=:blue, markersize = 15, marker = :diamond, label="Benchmark Jacobian + solver")

# ylims!(ax, 1.0, 5.0)

f

##

save("BM_MM_n62_to_n100.png", f)

##

# Intercept A: 3.619567894216926
# Slope A: 0.8196356820257618
# Error on Intercept A: 0.24105226536872967
# Error on Slope A: 0.12346154972613471
