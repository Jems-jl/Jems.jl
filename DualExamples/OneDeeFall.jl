#=
# OneDeeFall.jl

This notebook provides a simple example using dual numbers to compute the trajectory of a falling projectile in one dimension.
=#

import ForwardDiff: Dual
using GLMakie, LaTeXStrings, MathTeXEngine, Printf
basic_theme = Theme(fonts=(regular=texfont(:text), bold=texfont(:bold),
                           italic=texfont(:italic), bold_italic=texfont(:bolditalic)),
                    fontsize=20, size=(600, 450), px_per_unit=300, linewidth=4, markersize=5, scalefactor=4,
                    Axis=(xlabelsize=20, ylabelsize=20, titlesize=40, xgridvisible=false, ygridvisible=false,
                          spinewidth=2.5, xminorticksvisible=true, yminorticksvisible=true, xtickalign=1, ytickalign=1,
                          xminortickalign=1, yminortickalign=1, xticksize=14, xtickwidth=2.5, yticksize=14,
                          ytickwidth=2.5, xminorticksize=7, xminortickwidth=2.5, yminorticksize=7, yminortickwidth=2.5,
                          xticklabelsize=20, yticklabelsize=20, xticksmirrored=true, yticksmirrored=true),
                    Legend=(patchsize=(50, 10), framevisible=false, patchlabelgap=10, rowgap=5, labelsize=15))
set_theme!(basic_theme)
GLMakie.activate!(; scalefactor=4)
## convenience functions for working with dual numbers
to_value = x -> x.value
first_partial = x -> x.partials[1]
## define the system of ODEs
g = 1

dy_dt(y, vy) = vy
dvy_dt(y, vy) = -g


v0 = 0.0
t0 = 0.0
dt = dt_max = 0.3 # time step
dt_min = 0.01

y0s = LinRange(1.0, 10.0, 100)  # initial heights to test
Δy_target = 0.1
## define analytical solutions for comparison
y_analytic(t, y0) = y0 + v0 * t - (g / 2) * t^2
∂y_∂t_analytic(t, y0) = v0 - g * t
∂y_∂y₀_analytic(t, y0) = 1.0
T_analytic(y0) = sqrt(2 * y0 / g)
dT_dh_analytic(y, y0) = 1 / sqrt(2 * g * (y0 - y))

## setup arrays to store results for each initial angle
t_arrays = []
y_arrays = []
vy_arrays = []

function RK4(t, y, vy, dt)
    """
    Perform a single Runge-Kutta 4th order step for the system of ODEs defined by dx_dt, dvx_dt, dy_dt, dvy_dt.
    Returns the updated values of t, x, vx, y, vy after a time step
    """
    t_next = t + dt
    k1_y = dy_dt(y, vy)
    k1_vy = dvy_dt(y, vy)
    k2_y = dy_dt(y + 0.5 * dt * k1_y, vy + 0.5 * dt * k1_vy)
    k2_vy = dvy_dt(y + 0.5 * dt * k1_y, vy + 0.5 * dt * k1_vy)
    k3_y = dy_dt(y + 0.5 * dt * k2_y, vy + 0.5 * dt * k2_vy)
    k3_vy = dvy_dt(y + 0.5 * dt * k2_y, vy + 0.5 * dt * k2_vy)
    k4_y = dy_dt(y + dt * k3_y, vy + dt * k3_vy)
    k4_vy = dvy_dt(y + dt * k3_y, vy + dt * k3_vy)
    y_next = y + (dt / 6) * (k1_y + 2 * k2_y + 2 * k3_y + k4_y)
    vy_next = vy + (dt / 6) * (k1_vy + 2 * k2_vy + 2 * k3_vy + k4_vy)
    return t_next, y_next, vy_next
end

function forward_euler_step(t, y, vy, dt)
    """
    Perform a single forward Euler step for the system of ODEs defined by dx_dt, dvx_dt, dy_dt, dvy_dt.
    Returns the updated values of t, x, vx, y, vy after a time step
    """
    t_next = t + dt
    y_next = y + dy_dt(y, vy) * dt
    vy_next = vy + dvy_dt(y, vy) * dt
    return t_next, y_next, vy_next
end

for y0 in y0s
    y₀ = Dual(y0, 1.0)  # make y0 a dual number to compute derivatives with respect to it

    t_array = [Dual(t0, 0.0)]
    y_array = [Dual(y0, 1.0)]  # the dependent coordinates are also dual numbers to compute derivatives wrt y0
    vy_array = [Dual(v0, 0.0)]

    # do integration loop
    for i in 1:10000
        t = t_array[end]
        y = y_array[end]
        vy = vy_array[end]

        dt = min(Dual(dt_max, 0), max(abs(Δy_target / vy) * y.value / y₀, Dual(dt_min, 0)))  # adjust dt to keep the change in y small
        t_next, y_next, vy_next = RK4(t, y, vy, dt)
        push!(t_array, t_next)
        push!(y_array, y_next)
        push!(vy_array, vy_next)
        
        # Break if the object hits (or is in) the ground
        if y_next.value < 0.0
            break
        end
    end
    push!(t_arrays, t_array)
    push!(y_arrays, y_array)
    push!(vy_arrays, vy_array)
end
## plot the trajectories for each initial angle, along with the analytical solutions
f = Figure();
ax = Axis(f[1, 1], xlabel=L"t", ylabel=L"y")
for i in eachindex(y0s)
    scatter!(ax, to_value.(t_arrays[i]), to_value.(y_arrays[i]))
    lines!(ax, to_value.(t_arrays[i]), y_analytic.(to_value.(t_arrays[i]), Ref(y0s[i])), linewidth=3)
end
f
## define interpolation functions
function interpolate(x, x_array, y_array)  # interpolates y_array(x_array) at an input x
    for i in eachindex(x_array)[2:end]
        if x_array[i-1] < x <= x_array[i] || x_array[i] < x <= x_array[i-1]
            slope = (y_array[i] - y_array[i-1]) / (x_array[i] - x_array[i-1])
            return y_array[i-1] + (x - x_array[i-1]) * slope
        end
    end
    return Dual(NaN, NaN)
end

function slope(x, x_array, y_array)  # returns the slope of the linear interpolation of y_array(x_array) at an input x
    for i = 1:(length(x_array) - 1)
        if x_array[i] > x >= x_array[i+1] || x_array[i] < x <= x_array[i+1]
            return (y_array[i+1] - y_array[i]) / (x_array[i+1] - x_array[i])
        end
    end
    return Dual(NaN, NaN)
end

## find the t and y values at which the projectile hits the ground (y=0) for each initial angle
tends = []
yends = []
t_nearests = []
for i in 1:length(y0s)
    t_end = interpolate(0, y_arrays[i], t_arrays[i])
    y_end = interpolate(0, y_arrays[i], y_arrays[i])  # returns 0 by construction
    if abs(y_arrays[i][end].value) < abs(y_arrays[i][end-1].value)
        t_nearest = t_arrays[i][end]
    else
        t_nearest = t_arrays[i][end-1]
    end
    push!(t_nearests, t_nearest)
    push!(tends, t_end)
    push!(yends, y_end)
end

last_ders = []
second_to_last_ders = []
nearest_ders = []
true_ders = []
int_ders = []
for i in eachindex(y0s)
    push!(last_ders, first_partial.(t_arrays[i][end]))
    push!(second_to_last_ders, first_partial.(t_arrays[i][end-1]))
    push!(true_ders, dT_dh_analytic.(yends[i].value, y0s[i]))
    push!(int_ders, first_partial.(tends[i]))
    if abs(tends[i].value - t_arrays[i][end].value) < abs(tends[i].value - t_arrays[i][end-1].value)
        push!(nearest_ders, first_partial.(t_arrays[i][end]))
    else
        push!(nearest_ders, first_partial.(t_arrays[i][end-1]))
    end
end

##
##  display the time of fall as a function of y0, and the derivative of it with respect to y0
f = Figure();
ax = Axis(f[1, 1], ylabel=L"T", xticklabelsvisible=false)
lines!(ax, y0s, T_analytic.(y0s), color=:red, alpha=0.5, label="analytic")
scatter!(ax, y0s, to_value.(tends), color=:red,label="interpolated")
scatter!(ax, y0s, to_value.(t_nearests), color=:blue, label="nearest")
axislegend(ax, position=:rb)

ax2 = Axis(f[2, 1], xlabel=L"y_0")
# lines!(ax2, y0s, analytic.(t_end_analytic.(theta0s), theta0s), color=:blue, alpha=0.5,
    #    linestyle=:dash, label=L"\partial{X}/\partial{y_0}")
scatter!(ax2, y0s, nearest_ders, color=:blue, label="nearest")
lines!(ax2, y0s, true_ders, color=:red, alpha=0.5, linestyle=:dash, label=L"\mathrm{d}T/\mathrm{d}y_0")
scatter!(ax2, y0s, int_ders, color=:red, label="interpolated")
axislegend(ax2;)
linkxaxes!(ax, ax2)
f
save("DualExamples/T_fall_y0.png", f)

## choose the trajectory for y0=10
y_array = y_arrays[end]
t_array = t_arrays[end]
tend = tends[end]
yend = yends[end]

@inline write_to_string(x) = @sprintf("%0.2f", x)
dual_label(dualt, dualy) = L"\hat{t} = ({%$(write_to_string(dualt.value))}, {%$(write_to_string(dualt.partials[1]))}),
                             \hat{y} = ({%$(write_to_string(dualy.value))}, {%$(write_to_string(dualy.partials[1]))})"
int_dual_label(dualt, dualy) = L"\mathcal{\hat{T}} = ({%$(write_to_string(dualt.value))}, {%$(write_to_string(dualt.partials[1]))}),
                             \hat{\mathcal{Y}} = ({%$(write_to_string(dualy.value))}, {%$(write_to_string(dualy.partials[1]))})"

ts = LinRange(to_value(t_array[end-1]), to_value(t_array[end]), 100)
ys = interpolate.(Dual.(ts, 0), Ref(to_value.(t_array)), Ref(to_value.(y_array)))

## ...and plots its final points
f = Figure();
ax = Axis(f[1, 1], xlabel=L"t", ylabel=L"y")
scatter!(ax, to_value.(t_array)[end-3:end], to_value.(y_array)[end-3:end], color=:blue, label=L"y(t)", markersize=10)
lines!(ax, ts, to_value.(ys), label=L"\mathrm{interpolated} y(t)", color=:blue, linewidth=3, linestyle=:dash)
hlines!(ax, [0.0], color=:gray, linestyle=:dash, label=L"y=0")
scatter!(ax, to_value(tend), 0.0, color=:red, markersize=10)
ylims!(ax, -0.8, 3.8)
annotation!(ax, to_value(tend), 0.0, text=int_dual_label(tend, yend), align=(:right, :bottom), fontsize=15)
annotation!(ax, to_value(t_array[end]), to_value(y_array[end]), text=dual_label(t_array[end], y_array[end]), align=(:right, :top), fontsize=15)
annotation!(ax, to_value(t_array[end-1]), to_value(y_array[end-1]), text=dual_label(t_array[end-1], y_array[end-1]), align=(:left, :bottom), fontsize=15)
annotation!(ax, to_value(t_array[end-2]), to_value(y_array[end-2]), text=dual_label(t_array[end-2], y_array[end-2]), align=(:left, :bottom), fontsize=15)
annotation!(ax, to_value(t_array[end-3]), to_value(y_array[end-3]), text=dual_label(t_array[end-3], y_array[end-3]), align=(:left, :top), fontsize=15)
f
save("DualExamples/y_t_1D.png", f)