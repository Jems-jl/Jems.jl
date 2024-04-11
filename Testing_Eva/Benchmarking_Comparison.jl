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

# PLOT

data = x1.times
nbins = 30
minval, maxval = minimum(data), maximum(data)
intervals = LinRange(minval, maxval, nbins + 1)
steps = intervals[2] - intervals[1]
counts = [count(x -> intervals[i] <= x < intervals[i+1], data) for i in 1:nbins]
middle_values = [(intervals[i] + intervals[i+1]) / 2 for i in 1:length(intervals)-1]
middle_array = collect(middle_values)


barplot(
    middle_array, # x data
    counts, # y data
    color = "purple"

)


