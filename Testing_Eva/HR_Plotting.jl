
# This script is for plotting the values that are obtained for the different stellar models

# Load in the dataframes that are saved as csv files

using CSV
using DataFrames

# Kipp Rates

star_12M = CSV.read("history_12M_100R.csv", DataFrame)
star_10M = CSV.read("history_10M_100R.csv", DataFrame)
star_08M = CSV.read("history_08M_100R.csv", DataFrame)

L_surf_star_12 = log10.(star_12M.L_surf)
L_surf_star_10 = log10.(star_10M.L_surf)
L_surf_star_08 = log10.(star_08M.L_surf)

T_surf_star_12 = log10.(star_12M.T_surf)
T_surf_star_10 = log10.(star_10M.T_surf)
T_surf_star_08 = log10.(star_08M.T_surf)


# Toy Rates

star_10M_ToyRates = CSV.read("history_10M_100R_Toy.csv", DataFrame)
L_surf_star_10_Toy = log10.(star_10M_ToyRates.L_surf)
T_surf_star_10_Toy = log10.(star_10M_ToyRates.T_surf)

# Jina Rates

star_10M_JinaRates_V1 = CSV.read("history_10M_100R_Jina_Version_1.csv", DataFrame)
star_10M_JinaRates_V2 = CSV.read("history_10M_100R_Jina_Version_2.csv", DataFrame)

L_surf_star_10_Jina_1 = log10.(star_10M_JinaRates_V1.L_surf)
L_surf_star_10_Jina_2 = log10.(star_10M_JinaRates_V2.L_surf)

T_surf_star_10_Jina_1 = log10.(star_10M_JinaRates_V1.T_surf)
T_surf_star_10_Jina_2 = log10.(star_10M_JinaRates_V2.T_surf)

##

# Set up Makie

using CairoMakie, LaTeXStrings, MathTeXEngine
basic_theme = Theme(fonts=(regular=texfont(:text), bold=texfont(:bold),
                           italic=texfont(:italic), bold_italic=texfont(:bolditalic)),
                    fontsize=30, size=(1000, 750), linewidth=7,
                    Axis=(xlabelsize=40, ylabelsize=40, titlesize=40, xgridvisible=false, ygridvisible=false,
                          spinewidth=2.5, xminorticksvisible=true, yminorticksvisible=true, xtickalign=1, ytickalign=1,
                          xminortickalign=1, yminortickalign=1, xticksize=14, xtickwidth=2.5, yticksize=14,
                          ytickwidth=2.5, xminorticksize=7, xminortickwidth=2.5, yminorticksize=7, yminortickwidth=2.5,
                          xticklabelsize=35, yticklabelsize=35, xticksmirrored=true, yticksmirrored=true),
                    Legend=(patchsize=(70, 10), framevisible=false, patchlabelgap=20, rowgap=10))
set_theme!(basic_theme)

##

# Plot HR diagram

f = Figure();

ax = Axis(f[1, 1]; xlabel=L"\log_{10}(T_\mathrm{eff}/[K])", ylabel=L"\log_{10}(L/L_\odot)", xreversed=true)

# lines!(ax, T_surf_star_12, L_surf_star_12, color=:red, linewidth=3, label = L"1.2 M_\odot")
lines!(ax, T_surf_star_10, L_surf_star_10, color=:purple, linewidth=3, label = L"1.0 M_\odot\text{ Kippenhahn Rates}")
# lines!(ax, T_surf_star_08, L_surf_star_08, color=:blue, linewidth=3, label = L"0.8 M_\odot")

lines!(ax, T_surf_star_10_Toy, L_surf_star_10_Toy, color=:blue, linewidth=3, label = L"1.0 M_\odot \text{ Toy Rates}")

lines!(ax, T_surf_star_10_Jina_1, L_surf_star_10_Jina_1, color=:pink, linewidth=3, label = L"1.0 M_\odot \text{ Jina Rates version 1}")
# lines!(ax, T_surf_star_10_Jina_2, L_surf_star_10_Jina_2, color=:red, linewidth=3, label = L"1.0 M_\odot \text{ Jina Rates version 2}")


xlims!(ax, (4.3, 3.0)) 
# ylims!(ax, (ymin, ymax)) 


axislegend(position=:rt)

f

# save("comparison_Jina_Rates_Versions.png", f)
# save("comparison_Toy_Kipp_Jina_Rates.png", f)
# save("comparison_ToyRates_KippRates.png", f)
# save("comparison_masses_KippRates.png", f)



