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
import ForwardDiff.Dual

#defining differential equations
function interpolate(t_array, y_array, y0)
    for i in 1:length(y_array)-1
        if y_array[i] > y0 >= y_array[i+1]
            slope = (y_array[i+1] - y_array[i]) / (t_array[i+1] - t_array[i])
            return t_array[i] + (y0 - y_array[i]) / slope
        end
    end
    return NaN
end
function nearest_neighbour(t_array, y_array, y0)
    for i in 1:length(y_array)-1
        if y_array[i] > y0 >= y_array[i+1]
            if abs(y_array[i] - y0) < abs(y_array[i+1] - y0)
                return t_array[i]
            else
                return t_array[i+1]
            end
        end
    end
    return NaN
end

dy_dt(y,vy)  = vy
dvy_dt(y,vy) = -9.81 #m/s²
function free_fall(height::Float64)
    max_nb_steps = 100_000
    # defining initial conditions as dual numbers
    t0      = Dual(0.0,0.0)
    y0 = h  = Dual(height,1.0) # we keep track of derivatives with respect to height
    Delta_t = 0.1
    v0      = 0.0
    t_array = [t0]
    y_array = [y0]
    v_array = [v0]

    for i in 1:max_nb_steps
        # Forward Euler method
        t = t_array[end]
        y = y_array[end]
        v = v_array[end]
        t_next = t + Delta_t
        y_next = y + Delta_t * dy_dt(y, v)
        v_next = v + Delta_t * dvy_dt(y, v)
        push!(t_array, t_next)
        push!(y_array, y_next)
        push!(v_array, v_next)
        if y_next < 0.0
            break
        end
    end

    t_array
    y_array
    v_array
    slope = (y_array[end-1] - y_array[end-2]) / (t_array[end-1] - t_array[end-2])
    t_end = t_array[end]   - y_array[end]      /  slope
    T_interpo = interpolate(t_array, y_array,0)
    T_nearest = nearest_neighbour(t_array, y_array,0)
    return t_array, y_array, v_array, T_interpo, T_nearest
end

struct Simulation
    height
    t_array
    t_array_val
    t_array_partial
    y_array
    y_array_val
    y_array_partial
    v_array
    T_interpo
    T_nearest
end
function Simulation(height)
    (t_array, y_array, v_array, T_interpo, T_nearest) =  free_fall(height)
    t_array_val = [t.value for t in t_array]
    t_array_partial = [t.partials[1] for t in t_array]
    y_array_val = [y.value for y in y_array]
    y_array_partial = [y.partials[1] for y in y_array]
    Simulation(height, t_array, t_array_val, t_array_partial, y_array, y_array_val, y_array_partial, v_array, T_interpo, T_nearest)
end
struct Run
    simulations
    heights
    T_interpo
    T_interpo_val
    T_interpo_partial
    T_nearest
    T_nearest_val
    T_nearest_partial
    T_analytical_val
    T_analytical_partial
    diffs_interpo
    diffs_nearest
end

function Run(heights)
    simulations = [Simulation(height) for height in heights]
    T_interpo = [simulation.T_interpo for simulation in simulations]
    T_interpo_val = [T.value for T in T_interpo]
    T_interpo_partial = [T.partials[1] for T in T_interpo]
    T_nearest = [simulation.T_nearest for simulation in simulations]
    T_nearest_val = [T.value for T in T_nearest]
    T_nearest_partial = [T.partials[1] for T in T_nearest]
    T_analytical_val = [sqrt(2 * height / 9.81) for height in heights]
    T_analytical_partial = [1 / sqrt(2 * 9.81 * height) for height in heights]
    diffs_interpo = finite_diff(heights, T_interpo_val)
    diffs_nearest = finite_diff(heights, T_nearest_val)
    Run(simulations, heights, T_interpo, T_interpo_val, T_interpo_partial, T_nearest, T_nearest_val, T_nearest_partial, T_analytical_val, T_analytical_partial, diffs_interpo, diffs_nearest)
end

function plot(simulation::Simulation)
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel=L"t \, \mathrm{(s)}", ylabel=L"y \, \mathrm{(m)}")
    scatter!(ax, simulation.t_array_val, simulation.y_array_val; label="y(t)", color=:blue)
    scatter!(ax, simulation.T_interpo.value, 0; label="Interpolated", color=:red,markersize=20)
    scatter!(ax, simulation.T_nearest.value, 0; label="Nearest neighbour", color=:green, markersize=20,marker=:cross)
    axislegend()
    fig
end

function finite_diff(heights, Times)
    N = length(heights)
    diffs = zeros(N)
    for i in 1:N
        h = heights[i]
        t = Times[i]
        if i == 1
            Δh = heights[i+1] - heights[i]
            Δt = Times[i+1] - Times[i]
        else
            Δh = heights[i] - heights[i-1]
            Δt = Times[i] - Times[i-1]
        end
        diffs[i] = Δt / Δh
    end
    return diffs
end

##
heights = collect(0.3:0.15:15)
run = Run(heights)
##
fig = Figure(size = (1500, 1000))
ax = Axis(fig[1, 1]; ylabel=L"T(h)",xticklabelsvisible=false)
scatter!(ax, heights, run.T_interpo_val; label="Interpolated", color=:blue)
scatter!(ax, heights, run.T_nearest_val; label="Nearest neighbour", color=:red, marker=:cross, markersize=15)
lines!(ax, heights, run.T_analytical_val; label="Analytical", alpha = 0.8,linestyle=:dash, color=:orange)
axislegend(ax,position=:rb)
ax2 = Axis(fig[2, 1];xticklabelsvisible=true,ylabel=L"\partial T / \partial h",xlabel=L"h \, \mathrm{(m)}")
linkxaxes!(ax, ax2)
#scatter!(ax2, heights, run.diffs_interpo; label="Interpolated", color=:blue)
scatter!(ax2, heights, run.T_interpo_partial; label="Partial from interpolated dual", color=:blue)
#scatter!(ax2, heights, run.diffs_nearest; label=L"Nearest neighbour $(T_{i+1}-T_i)/(h_{i+1}-h_i)$", color=:red, marker=:cross, markersize=15)
scatter!(ax2, heights, run.T_nearest_partial; label="Partial from nearest neighbour dual", color=:red, marker=:cross, markersize=15)
#lines!(ax2, heights, run.diffs_nearest; color=:red,alpha=0.1)
lines!(ax2, heights, run.T_analytical_partial; alpha=0.8,linestyle=:dash,label="Analytical partial", color=:orange)
#scatter!(ax2, heights, run.T_interpo_partial; label="Interpolated partial", color=:orange)
axislegend(ax2)
#leg = Legend(fig[0,1], ax2, position=:lt,fontsize=20,orientation=:horizontal,tellwidth =false)
fig
##
savepath = "Figures/partial_plot_pm.png"; save(savepath, fig; px_per_unit = 3); @show savepath
##