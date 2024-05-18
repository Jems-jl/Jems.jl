using Interpolations
using CairoMakie
import ForwardDiff.Dual
##
# Lower and higher bound of interval
a = 1.0
b = 10.0
# Interval definition
x = a:1.0:b
# This can be any sort of array data, as long as
# length(x) == length(y)
y = @. cos(x^2 / 9.0) # Function application by broadcasting
# Interpolations
itp_linear = linear_interpolation(x, y)
itp_cubic = cubic_spline_interpolation(x, y)
itp_cubic = interpolate(x,y,FritschCarlsonMonotonicInterpolation())
# Interpolation functions
f_linear(x) = itp_linear(x)
f_cubic(x) = itp_cubic(x)
function finite_diff(x,y)
    N = length(x)
    h = x[2] - x[1]
    D = zeros(N)
    for i in 1:N
        if i == 1
            D[i] = (y[i+1] - y[i])/h
        elseif i == N
            D[i] = (y[i] - y[i-1])/h
        else
            D[i] = (y[i+1] - y[i-1])/(2h)
        end
    end
    return D
end

# Plots
fig = Figure(resolution = (1200, 800))
ax = Axis(fig[1, 1], xlabel = "x", ylabel = "f(x)")

x_new = collect(a:0.01:b) # smoother interval, necessary for cubic spline
y_linear = f_linear.(x_new)
y_cubic = f_cubic.(x_new)
diff_linear = finite_diff(x_new, y_linear)
diff_cubic = finite_diff(x_new, y_cubic)
diff2_cubic = finite_diff(x_new, diff_cubic)
scatter!(ax, x, y, markersize = 10, color = :black, label = "Data points")
lines!(ax, x_new, y_linear, color = :red, linewidth = 1, label = "Linear interpolation")
lines!(ax, x_new, diff_linear, markersize = 3, color = :red, label = "Linear derivative")

lines!(ax, x_new, y_cubic, color = :blue, linewidth = 1, linestyle=:dash, label = "Cubic interpolation")
lines!(ax, x_new, diff_cubic, markersize = 3,  color = :blue, label =  "Cubic Derivative")
lines!(ax, x_new, diff2_cubic, markersize = 3, color = :blue,linestyle=:dot, label = "Cubic Second Derivative")
axislegend(ax)
fig