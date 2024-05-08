using BenchmarkTools
using LinearAlgebra
using Statistics

# reload the values from the benchmarks

x_start = 4
x_end   = 300
n_data  = 25

x_flt = 10 .^(range(log10(x_start), stop = log10(x_end), length = n_data))      # x_flt = float values of normal x-values
x_rnd = round.(Int, x_flt)                                                      # x_rnd = rounded x-values in normal scale
x_log = log10.(x_rnd)                                                           # x_log = x-values in log scale


folder_path = "/Users/evakuipers/Downloads/All-basic-data_2024-05-07_1523/BM_LDIV_SA_output"
benchmark_files = readdir(folder_path)
json_files = filter(x -> endswith(x, ".json"), benchmark_files)


function custom_sort(filename1::String, filename2::String)
    number1 = parse(Int, match(r"\d+", filename1).match)
    number2 = parse(Int, match(r"\d+", filename2).match)
    return number1 < number2
end

sorted_files = sort(json_files, lt=custom_sort)


m_list = []
s_list = []
q10_list = []
q25_list = []
q75_list = []
q90_list = []

for file in sorted_files

    bm = BenchmarkTools.load(joinpath(folder_path, file))[1]
    data = bm.times
    m = median(data); push!(m_list, m)
    s = std(bm).time; push!(s_list, s)

    perc_10 = quantile(data, 0.10); push!(q10_list, perc_10)
    perc_25 = quantile(data, 0.25); push!(q25_list, perc_25)
    perc_75 = quantile(data, 0.75); push!(q75_list, perc_75)
    perc_90 = quantile(data, 0.90); push!(q90_list, perc_90)
    
end


y_err_top = m_list .- q25_list 
y_err_bot = q75_list .- m_list

x_values = x_rnd
y_values = m_list

y_errors_top = y_err_top ./ m_list
y_errors_bot = y_err_bot ./ m_list

y_errors_s = s_list ./ m_list


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


# function weighted_least_squares_fit(xdata, ydata, yerrors)
#     @assert length(xdata) == length(ydata) == length(yerrors)
#     n = length(xdata)
#     X = hcat(ones(n), xdata) 
#     W = Diagonal(1 ./ (yerrors .^ 2))
#     β = W * X \ (W * ydata)
#     cov_matrix = inv(X' * W * X)
#     return β, cov_matrix
# end

# α, cov_matrix_A = weighted_least_squares_fit(log10.(x_values), log10.(y_values), y_errors_s)

# function fit_errors(cov_matrix)
#     errors = sqrt.(diag(cov_matrix))
#     return errors
# end

# println("Intercept A: ", α[1])
# println("Slope A: ", α[2])
# println("Error on Intercept A: ", fit_errors(cov_matrix)[1])
# println("Error on Slope A: ", fit_errors(cov_matrix)[2])

# fitted_line_A(x) = α[1] + α[2] * x


function least_squares_fit(xdata, ydata)
    @assert length(xdata) == length(ydata)
    n = length(xdata)
    X = hcat(ones(n), xdata) 
    β = X \ ydata
    cov_matrix = inv(X' * X) 
    return β, cov_matrix
end

α, cov_matrix_A = least_squares_fit(log10.(x_values), log10.(y_values))

function fit_errors(cov_matrix)
    errors = sqrt.(diag(cov_matrix))
    return errors
end

println("*** Values for LDIV - using StaticArrays ***")
println("Intercept A: ", α[1])
println("Slope A: ", α[2])
println("Error on Intercept A: ", fit_errors(cov_matrix_A)[1])
println("Error on Slope A: ", fit_errors(cov_matrix_A)[2])

fitted_line_A(x) = α[1] + α[2] * x


f = Figure();
ax = Axis(f[1, 1]; xlabel=L"\log_{10}(\mathrm{n})", ylabel=L"\log_{10}(\mathrm{t/[ns]})")

lines!(log10.(x_values), x -> fitted_line_A(x), color=:grey, linewidth=2, linestyle=:solid, label="Fit")
errorbars!(log10.(x_values), log10.(y_values), y_errors_top, y_errors_bot, color=:blue, linewidth=2, whiskerwidth=10)
scatter!(log10.(x_values), log10.(y_values), color=:blue, markersize = 15, marker = :diamond, label="Benchmark Jacobian + solver")

f

##

# save("BM_MM_n62_to_n100.png", f)

