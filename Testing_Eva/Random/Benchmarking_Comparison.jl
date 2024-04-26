using BenchmarkTools
using Makie
using CairoMakie


##

# RUN1 --> [:H1, :He4, :C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
# RUN2 --> [:H1, :D2, :He4, :C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
# RUN3 --> [:H1, :D2, :He3, :He4, :C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
# RUN4 --> [:H1, :D2, :He3, :He4, :Li7, :C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
# RUN5 --> [:H1, :D2, :He3, :He4, :Li7, :Be8, :C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])


# Load Files

x1_file = BenchmarkTools.load("Kipp_RUN1_a.json"); x1 = x1_file[1]
x2_file = BenchmarkTools.load("Kipp_RUN2_a.json"); x2 = x2_file[1]
x3_file = BenchmarkTools.load("Kipp_RUN3_a.json"); x3 = x3_file[1]
x4_file = BenchmarkTools.load("Kipp_RUN4_a.json"); x4 = x4_file[1]
x5_file = BenchmarkTools.load("Kipp_RUN5_a.json"); x5 = x5_file[1]

y1_file = BenchmarkTools.load("Kipp_RUN1_b.json"); y1 = y1_file[1]
y2_file = BenchmarkTools.load("Kipp_RUN2_b.json"); y2 = y2_file[1]
y3_file = BenchmarkTools.load("Kipp_RUN3_b.json"); y3 = y3_file[1]
y4_file = BenchmarkTools.load("Kipp_RUN4_b.json"); y4 = y4_file[1]
y5_file = BenchmarkTools.load("Kipp_RUN5_b.json"); y5 = y5_file[1]


Kipp_file_a = BenchmarkTools.load("Kipp_RUN1_a.json"); Kipp_a = Kipp_file_a[1]
Jina_file_a = BenchmarkTools.load("Jina_RUN1_a.json"); Jina_a = Jina_file_a[1]

Kipp_file_b = BenchmarkTools.load("Kipp_RUN1_b.json"); Kipp_b = Kipp_file_b[1]
Jina_file_b = BenchmarkTools.load("Jina_RUN1_b.json"); Jina_b = Jina_file_b[1]



# Start comparing the benchmarks

##
BenchmarkTools.judge(median(x1), median(x2))
##
BenchmarkTools.judge(median(x2), median(x3))
##
BenchmarkTools.judge(median(x3), median(x4))
##
BenchmarkTools.judge(median(x4), median(x5))
##


##
BenchmarkTools.judge(median(y1), median(y2))
##
BenchmarkTools.judge(median(y2), median(y3))
##
BenchmarkTools.judge(median(y3), median(y4))
##
BenchmarkTools.judge(median(y4), median(y5))
##


# Compare Jina and Kipp

BenchmarkTools.judge(median(Kipp_a), median(Jina_a))
##
BenchmarkTools.judge(median(Kipp_b), median(Jina_b))
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

# PLOT

f = Figure()

ax1 = Axis(f[1, 1])

data = x1.times
nbins = 30
minval, maxval = minimum(data), maximum(data)
intervals = LinRange(minval, maxval, nbins + 1)
steps = intervals[2] - intervals[1]
counts = [count(x -> intervals[i] <= x < intervals[i+1], data) for i in 1:nbins]
middle_values = [(intervals[i] + intervals[i+1]) / 2 for i in 1:length(intervals)-1]
middle_array = collect(middle_values)

barplot!(
    middle_array,       # x data
    counts,             # y data
    color = :orchid1      # color

)
# vlines!([median(x1).time,
#          median(x1).time - std(x1).time,
#          median(x1).time + std(x1).time],
#          linewidth = 3,
#          color = :blue)


f
