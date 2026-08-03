import ForwardDiff: Dual
using NLsolve
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
##
to_value = x -> x.value
first_partial = x -> x.partials[1]
second_partial = x -> x.partials[2]
third_partial = x -> x.partials[3]
## define the system of ODEs for projectile motion
g = 1

dx_dt(x, vx, y, vy) = vx
dvx_dt(x, vx, y, vy) = 0
dy_dt(x, vx, y, vy) = vy
dvy_dt(x, vx, y, vy) = -g

v0 = 1.0
y0 = 0.0
x0 = 0.0

t0 = 0.0
theta0s = LinRange(pi/8, 3*pi/8, 100)
dt = 0.08 # time step
##
function x_range_analytic(theta0)
    return x0 + v0^2 * cos(theta0) * sin(theta0) / g * (1 + sqrt(1 + 2*g*(y0)/v0^2/sin(theta0)^2))
end

function y_path_analytic(x, theta0)
    y0 + (x - x0) * tan(theta0) - (g / (2 * v0^2 * cos(theta0)^2)) * (x - x0)^2
end
y_analytic(t) = y0 + v0 * sin(theta0.value) * t - (g / 2) * t^2
x_analytic(t) = x0 + v0 * cos(theta0.value) * t
∂y_∂t_analytic(t) = v0 * sin(theta0.value) - g * t
∂x_∂t_analytic(t) = v0 * cos(theta0.value)
∂y_∂v₀_analytic(t) = sin(theta0.value) * t
∂y_∂θ₀_analytic(t) = v0 * cos(theta0.value) * t
∂y_∂₀_analytic(t) = 1.0
∂x_∂v₀_analytic(t) = cos(theta0.value) * t
∂x_∂θ₀_analytic(t, theta0) = -v0 * sin(theta0) * t
dx_dθ₀_analytic(theta0) = 2*cos(2*theta0)
t_end_analytic(theta0) = v0 * sin(theta0) / g * (1 + sqrt(1+ 2 * g * y0 / (v0^2 * sin(theta0)^2)))

function x_path_analytic(y)
    if y > y0 + v0^2 * sin(theta0.value)^2 / (2*g)
        return NaN
    end
    return x0 + v0^2 * cos(theta0.value) * sin(theta0.value) / g +
           v0^2 * cos(theta0.value)^2 / g *
           sqrt(tan(theta0.value)^2 + 2*g*(y0 - y)/v0^2/cos(theta0.value)^2)
end
function dx_dtheta_path_analytic(y, theta0)
    if y > y0 + v0^2 * sin(theta0)^2 / (2*g)
        return NaN
    end
    sq = sqrt(1 + 2*g*(y0 - y)/v0^2/sin(theta0)^2)
    return v0^2 * cos(2*theta0) / g * (1+sq) - (y0 - y) / sq * 2 * cot(theta0)^2
end

##
t_arrays = []
y_arrays = []
x_arrays = []
vx_arrays = []
vy_arrays = []

function RK4(t, x, vx, y, vy, dt)
    t_next = t + dt
    k1_x = dx_dt(x, vx, y, vy)
    k1_vx = dvx_dt(x, vx, y, vy)
    k1_y = dy_dt(x, vx, y, vy)
    k1_vy = dvy_dt(x, vx, y, vy)
    k2_x = dx_dt(x + 0.5 * dt * k1_x, vx + 0.5 * dt * k1_vx, y + 0.5 * dt * k1_y, vy + 0.5 * dt * k1_vy)
    k2_vx = dvx_dt(x + 0.5 * dt * k1_x, vx + 0.5 * dt * k1_vx, y + 0.5 * dt * k1_y, vy + 0.5 * dt * k1_vy)
    k2_y = dy_dt(x + 0.5 * dt * k1_x, vx + 0.5 * dt * k1_vx, y + 0.5 * dt * k1_y, vy + 0.5 * dt * k1_vy)
    k2_vy = dvy_dt(x + 0.5 * dt * k1_x, vx + 0.5 * dt * k1_vx, y + 0.5 * dt * k1_y, vy + 0.5 * dt * k1_vy)
    k3_x = dx_dt(x + 0.5 * dt * k2_x, vx + 0.5 * dt * k2_vx, y + 0.5 * dt * k2_y, vy + 0.5 * dt * k2_vy)
    k3_vx = dvx_dt(x + 0.5 * dt * k2_x, vx + 0.5 * dt * k2_vx, y + 0.5 * dt * k2_y, vy + 0.5 * dt * k2_vy)
    k3_y = dy_dt(x + 0.5 * dt * k2_x, vx + 0.5 * dt * k2_vx, y + 0.5 * dt * k2_y, vy + 0.5 * dt * k2_vy)
    k3_vy = dvy_dt(x + 0.5 * dt * k2_x, vx + 0.5 * dt * k2_vx, y + 0.5 * dt * k2_y, vy + 0.5 * dt * k2_vy)
    k4_x = dx_dt(x + dt * k3_x, vx + dt * k3_vx, y + dt * k3_y, vy + dt * k3_vy)
    k4_vx = dvx_dt(x + dt * k3_x, vx + dt * k3_vx, y + dt * k3_y, vy + dt * k3_vy)
    k4_y = dy_dt(x + dt * k3_x, vx + dt * k3_vx, y + dt * k3_y, vy + dt * k3_vy)
    k4_vy = dvy_dt(x + dt * k3_x, vx + dt * k3_vx, y + dt * k3_y, vy + dt * k3_vy)
    x_next = x + (dt / 6) * (k1_x + 2 * k2_x + 2 * k3_x + k4_x)
    vx_next = vx + (dt / 6) * (k1_vx + 2 * k2_vx + 2 * k3_vx + k4_vx)
    y_next = y + (dt / 6) * (k1_y + 2 * k2_y + 2 * k3_y + k4_y)
    vy_next = vy + (dt / 6) * (k1_vy + 2 * k2_vy + 2 * k3_vy + k4_vy)
    return t_next, x_next, vx_next, y_next, vy_next
end

function forward_euler_step(t, x, vx, y, vy, dt)
    t_next = t + dt
    x_next = x + dx_dt(x, vx, y, vy) * dt
    vx_next = vx + dvx_dt(x, vx, y, vy) * dt
    y_next = y + dy_dt(x, vx, y, vy) * dt
    vy_next = vy + dvy_dt(x, vx, y, vy) * dt
    return t_next, x_next, vx_next, y_next, vy_next
end

for theta0 in theta0s
    θ₀ = Dual(theta0, 1.0)

    t_array = [t0]
    y_array = [Dual(y0, 0.0)]
    x_array = [Dual(x0, 0.0)]
    vx_array = [v0 * cos(θ₀)]
    vy_array = [v0 * sin(θ₀)]

    ## do integration loop

    for i in 1:10000
        # Forward Euler method
        t = t_array[end]
        x = x_array[end]
        y = y_array[end]
        vx = vx_array[end]
        vy = vy_array[end]

        t_next, x_next, vx_next, y_next, vy_next = RK4(t, x, vx, y, vy, dt)
        # Store the results
        push!(t_array, t_next)
        push!(x_array, x_next)
        push!(y_array, y_next)
        push!(vx_array, vx_next)
        push!(vy_array, vy_next)
        
        # Break if the object hits the ground
        if y_next.value < 0.0
            break
        end
    end
    push!(t_arrays, t_array)
    push!(x_arrays, x_array)
    push!(y_arrays, y_array)
    push!(vx_arrays, vx_array)
    push!(vy_arrays, vy_array)
end
##
f = Figure();
ax = Axis(f[1, 1])
for i in eachindex(theta0s)
    scatter!(ax, to_value.(x_arrays[i]), to_value.(y_arrays[i]))
    lines!(ax, to_value.(x_arrays[i]), y_path_analytic.(to_value.(x_arrays[i]), theta0s[i]), linewidth=3)
end
f
## define interpolation and analytical functions for comparison

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
function find_t_end(t_array, y_array)
    return interpolate(0, y_array, t_array)
end
function find_x_end(x_array, y_array)
    return interpolate(0, y_array, x_array)
end

##
xends = []
yends = []
x_nearests = []
for i in 1:length(theta0s)
    x_end = interpolate(0, y_arrays[i], x_arrays[i])
    y_end = interpolate(0, y_arrays[i], y_arrays[i])
    if abs(x_end.value - x_arrays[i][end].value) < abs(x_end.value - x_arrays[i][end-1].value)
        x_nearest = x_arrays[i][end]
    else
        x_nearest = x_arrays[i][end-1]
    end
    push!(x_nearests, x_nearest)
    push!(xends, x_end)
    push!(yends, y_end)
end

last_ders = []
second_to_last_ders = []
nearest_ders = []
true_ders = []
int_ders = []
for i in eachindex(theta0s)
    push!(last_ders, first_partial.(x_arrays[i][end]))
    push!(second_to_last_ders, first_partial.(x_arrays[i][end-1]))
    push!(true_ders, dx_dtheta_path_analytic.(0.0, theta0s[i]))
    push!(int_ders, first_partial.(xends[i]))
    if abs(xends[i].value - x_arrays[i][end].value) < abs(xends[i].value - x_arrays[i][end-1].value)
        push!(nearest_ders, first_partial.(x_arrays[i][end]))
    else
        push!(nearest_ders, first_partial.(x_arrays[i][end-1]))
    end
end

##
##  display the x_range as a function of theta0, and the derivative of x_range with respect to theta0
f = Figure();
ax = Axis(f[1, 1], xticks=LinearTicks(6), ylabel=L"x", xticklabelsvisible=false)
lines!(ax, rad2deg.(theta0s), x_range_analytic.(theta0s), color=:red, alpha=0.5, label="analytic")
scatter!(ax, rad2deg.(theta0s), to_value.(xends), color=:red,label="interpolated")
# scatter!(ax, rad2deg.(theta0s), to_value.([x_arrays[i][end] for i = 1:length(theta0s)]), label="last step")
# scatter!(ax, rad2deg.(theta0s), to_value.([x_arrays[i][end-1] for i = 1:length(theta0s)]), label="second to last step")
scatter!(ax, rad2deg.(theta0s), to_value.(x_nearests), color=:blue, label="nearest")
axislegend(ax, position=:cb)

ax2 = Axis(f[2, 1], xticks=LinearTicks(6), xlabel=L"\theta_0")
# scatter!(ax, rad2deg.(theta0s), last_ders, color=:blue, label="last step")
# scatter!(ax, rad2deg.(theta0s), second_to_last_ders, color=:orange, label="second to last step")
lines!(ax2, rad2deg.(theta0s), ∂x_∂θ₀_analytic.(t_end_analytic.(theta0s), theta0s), color=:blue, alpha=0.5,
       linestyle=:dash, label=L"\partial{X}/\partial{\theta_0}")
scatter!(ax2, rad2deg.(theta0s), nearest_ders, color=:blue, label="nearest")
lines!(ax2, rad2deg.(theta0s), true_ders, color=:red, alpha=0.5, linestyle=:dash, label=L"\mathrm{d}X/\mathrm{d}\theta_0")
scatter!(ax2, rad2deg.(theta0s), int_ders, color=:red, label="interpolated")
axislegend(ax2;)
linkxaxes!(ax, ax2)
f
save("developing/range_theta_0.png", f)

##

x_array = x_arrays[div(length(theta0s),2)]
y_array = y_arrays[div(length(theta0s),2)]
t_array = t_arrays[div(length(theta0s),2)]
xend = xends[div(length(theta0s),2)]
yend = yends[div(length(theta0s),2)]

@inline write_to_string(x) = @sprintf("%0.2f", x)
dual_label(dualx, dualy) = L"\hat{x} = ({%$(write_to_string(dualx.value))}, {%$(write_to_string(dualx.partials[1]))}),
                             \hat{y} = ({%$(write_to_string(dualy.value))}, {%$(write_to_string(dualy.partials[1]))})"
int_dual_label(dualx, dualy) = L"\mathcal{\hat{X}} = ({%$(write_to_string(dualx.value))}, {%$(write_to_string(dualx.partials[1]))}),
                             \hat{\mathcal{Y}} = ({%$(write_to_string(dualy.value))}, {%$(write_to_string(dualy.partials[1]))})"

xs = LinRange(to_value(x_array[end-1]), to_value(x_array[end]), 100)
ys = interpolate.(Dual.(xs, 0), Ref(to_value.(x_array)), Ref(to_value.(y_array)))

##
f = Figure();
ax = Axis(f[1, 1], xlabel=L"x", ylabel=L"y")
scatter!(ax, to_value.(x_array)[end-3:end], to_value.(y_array)[end-3:end], color=:blue, label=L"y(x)", markersize=10)
lines!(ax, xs, to_value.(ys), label=L"\mathrm{interpolated} y(x)", color=:blue, linewidth=3, linestyle=:dash)
hlines!(ax, [0.0], color=:gray, linestyle=:dash, label=L"y=0")
scatter!(ax, to_value(xend), 0.0, color=:red, markersize=10)
annotation!(ax, to_value(xend), 0.0, text=int_dual_label(xend, yend), align=(:right, :top), fontsize=15)
annotation!(ax, to_value(x_array[end]), to_value(y_array[end]), text=dual_label(x_array[end], y_array[end]), align=(:right, :bottom), fontsize=15)
annotation!(ax, to_value(x_array[end-1]), to_value(y_array[end-1]), text=dual_label(x_array[end-1], y_array[end-1]), align=(:left, :bottom), fontsize=15)
annotation!(ax, to_value(x_array[end-2]), to_value(y_array[end-2]), text=dual_label(x_array[end-2], y_array[end-2]), align=(:left, :bottom), fontsize=15)
annotation!(ax, to_value(x_array[end-3]), to_value(y_array[end-3]), text=dual_label(x_array[end-3], y_array[end-3]), align=(:left, :top), fontsize=15)


f