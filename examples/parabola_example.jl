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

import ForwardDiff.Dual

#defining differential equations
dy_dt(y,vy)  = vy
dvy_dt(y,vy) = -9.81 #m/s²

max_nb_steps = 100_000

# defining initial conditions as dual numbers
t0      = Dual(0.0  ,0.0)

y0 = h  = Dual(10.0 ,1.0) # we keep track of derivatives with respect to height
Delta_t = h/1000 #Dual(0.1,0.0)
v0      = Dual(0.0  ,0.0) 
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
t_end = t_array[end-1] - y_array[end-1] /   slope
t_array[end]
function interpolate(t_array, y_array, y)
    for i in 1:length(y_array)-1
        if y_array[i] > y >= y_array[i+1]
            slope = (y_array[i+1] - y_array[i]) / (t_array[i+1] - t_array[i])
            return t_array[i] + (y - y_array[i]) / slope
        end
    end
    return NaN
end
interpolate(t_array, y_array, 0)
analytical_dt_dh(y,h) = 1 / sqrt(2 * 9.81 * (h-y))
analytical_v(y,h) = sqrt(2 * 9.81 * (h-y))
analytical_dt_dy(y,h) = - 1 / sqrt(2*9.81*(h-y))
analytical_t(y,h) = sqrt(2 * (h-y) / 9.81)
#y_interps = LinRange(y_array[1].value, y_array[end].value, 900)
y_interps = LinRange(y_array[3].value, y_array[end].value, 2000)
t_interps = interpolate.(Ref(t_array), Ref(y_array), y_interps)
##
f = Figure();
ax = Axis(f[1,1])
lines!(ax, y_interps,(d->d.value).(t_interps))
scatter!(ax, (d->d.value).(y_array),(d->d.value).(t_array), color = :red)
ax.xlabel = "t"
ax.ylabel = "y(t)"
f
##

f = Figure();
ax = Axis(f[1,1], xlabel= L"Particle height $y$", ylabel= L"\partial t / \partial h")
lines!(ax, y_interps,(d -> d.partials[1]).(t_interps), color = :blue, size = 0.1, label="Interpolation")
lines!(ax, y_interps,analytical_dt_dh.(y_interps, 10), color = :red, linewidth = 2, linestyle = :dash, label="Analytical")
axislegend()


f
##
fig = Figure(size=(500, 600))
ax1 = Axis(fig[1,1],xreversed = true,ylabel=L"Time $t$ [s]")
ax4 = Axis(fig[4,1], xreversed = true,ylabel=L"$\partial t / \partial h$ [s/m]", xlabel=L"Particle height $y$")
sample = scatter!(ax1, (d->d.value).(y_array),(d->d.value).(t_array), color = :red,label=L"$t_i(y_i)$")
inter = lines!(ax1, (d->d.value).(y_array),(d->d.value).(t_array), color = :blue,label=L" $t(y)$" )
analyt = lines!(ax1, y_interps, analytical_t.(y_interps, 10), color = :lime, linewidth = 2, linestyle = :dash, label=L"$t(y)$")
lines!(ax4, y_interps,(d -> d.partials[1]).(t_interps), color = :blue, size = 0.1, label=L"$\partial t / \partial h$")
lines!(ax4, y_interps,analytical_dt_dh.(y_interps, 10), color = :lime, linewidth = 2, linestyle = :dash, label=L"$\partial t / \partial h$")
ax1.xticks = (d->d.value).(y_array[2:end])
ax1.ygridvisible = false
ax4.xticks = (d->d.value).(y_array[2:end])
hidexdecorations!(ax1,grid=false);hidexdecorations!(ax4,grid=false)
ax4_dummy = Axis(fig[4, 1],xreversed=true,xticks = 10:-1:0, xlabel=L"Particle height $y$ [m]",yticks = [0,10])
hidespines!(ax4_dummy);hidexdecorations!(ax4_dummy)
ax4.ygridvisible = false
ax4_dummy.ygridvisible = false
ax4_dummy.yticksvisible = false
ax4_dummy.yticklabelsvisible = false  # = [0,10]
ax4_dummy.xticklabelsvisible = true
ax4_dummy.xticksvisible = true
ax4_dummy.xlabelvisible = true
linkyaxes!(ax4, ax4_dummy);linkxaxes!(ax4,  ax4_dummy);linkxaxes!(ax1, ax4)

ax2 = Axis(fig[3, 1], ylabel=L"$\partial t / \partial y$ [s/m]",xreversed=true,ygridvisible=false,xticksvisible=false,xticklabelsvisible=false,xticks = (d->d.value).(y_array[2:end]))
linkxaxes!(ax2, ax1)


interp = lines!(ax2, y_interps, (d->d.partials[1]).(interpolate.(Ref(t_array), Ref((d->d.value).(y_array)), (v->Dual(v,1.0)).(y_interps))), color = :blue, size = 0.1, label=L"$\partial t/\partial y$ ")
analyty = lines!(ax2, y_interps,analytical_dt_dy.(y_interps, 10), color = :lime, linewidth = 2, linestyle = :dash, label=L"$\partial t / \partial y$")
axislegend(ax2, titlefont=:regular, orientation=:horizontal,position=:rb,[interp,analyty],[L"$\partial t/\partial y$",L"$\partial t/\partial y$"], "Dual information from interpolation")

axislegend(ax1, orientation=:horizontal,position=:rb)

axislegend(ax4, orientation=:horizontal,position=:rt)
leg = Legend(fig[0,1],[sample, inter, analyt], ["Sampled points", "Linear interpolation", "Analytical result"],orientation=:horizontal,tellheight = true, tellwidth=true)
ax3 = Axis(fig[2,1], xreversed=true,ygridvisible=false,xticksvisible=false,xticklabelsvisible=false,xticks = (d->d.value).(y_array[2:end])); linkxaxes!(ax3, ax4)
dot = scatter!(ax3, (d->d.value).(y_array) , (d->d.partials[1]).(y_array), color = :red,label=L"$\partial y_i / \partial h$")
cross= scatter!(ax3, (d->d.value).(y_array),(d->d.partials[1]).(t_array), color = :red,label=L"$\partial t_i / \partial h$",marker = :x)
l = axislegend(ax3, titlefont=:regular, orientation = :horizontal,[dot,cross],[L"$\partial y_i / \partial h$", L"$\partial t_i / \partial h$ [s/m]"], "Dual information from simulation", position= :rc)
#plot vertical line
vlines!(ax1, [0], color = :black, linestyle = :dash, linewidth = 2)
vlines!(ax2, [0], color = :black, linestyle = :dash, linewidth = 2)
vlines!(ax3, [0], color = :black, linestyle = :dash, linewidth = 2)
vlines!(ax4, [0], color = :black, linestyle = :dash, linewidth = 2)
fig
##save Figure
save("parabola_example_master.png", fig,px_per_unit=11)
##
import ForwardDiff.Dual

# Define differential equations
dy_dt(y, vy) = vy
dvy_dt(y, vy) = -9.81 # m/s²

# Maximum number of steps
max_nb_steps = 100_00000

# Define initial conditions as dual numbers
t0 = Dual(0.0, 0.0)
y0 = h = Dual(10.0, 1.0) # Initial height
v0 = Dual(0.0, 0.0)
t_array = [t0]
y_array = [y0]
v_array = [v0]

# Define parameters for adaptive time stepping
target_dy = 0.01  # Target change in y per time step
min_dt = 0.0000001   # Minimum time step size

for i in 1:max_nb_steps
    # Forward Euler method with adaptive time step
    t = t_array[end]
    y = y_array[end]
    v = v_array[end]

    
    # Estimate the next time step size based on the rate of change of y
    dy = abs(dy_dt(y, v))
    Delta_t = min(max(target_dy / dy, min_dt), 0.1)  # Adjusted time step
    
    # Update variables using the Euler method
    t_next = t + Delta_t
    y_next = y + Delta_t * dy_dt(y, v)
    v_next = v + Delta_t * dvy_dt(y, v)
    
    # Store the results
    push!(t_array, t_next)
    push!(y_array, y_next)
    push!(v_array, v_next)
    
    # Break if the object hits the ground
    if y_next.value < 0.0
        break
    end
end



t_array
y_array
v_array
##
import ForwardDiff.Dual

#defining differential equations
k = 1
dx_dt(x,v)  = v
dv_dt(x,v) = -Dual(k,1)*x

max_nb_steps = 100_000

# defining initial conditions as dual numbers
t0      = Dual(0.0  ,0.0)
Delta_t = Dual(0.001,0.0)
x0      = Dual(10.0 ,1.0) # we keep track of derivatives with respect to height
v0      = Dual(0.0  ,0.0) 
t_array = [t0]
x_array = [x0]
v_array = [v0]

for i in 1:max_nb_steps
    # Forward Euler method
    t = t_array[end]
    x = x_array[end]
    v = v_array[end]
    t_next = t + Delta_t
    x_next = x + Delta_t * dx_dt(x, v)
    v_next = v + Delta_t * dv_dt(x, v)
    push!(t_array, t_next)
    push!(x_array, x_next)
    push!(v_array, v_next)
    if x_next < 0.0
        break
    end
end

t_array
x_array
v_array

slope = (x_array[end-1] - x_array[end-2]) / (t_array[end-1] - t_array[end-2])
t_end = t_array[end]   - x_array[end]    /  slope
t_end = t_array[end-1] - x_array[end-1] /   slope
t_array[end]

k = Dual(1.0, 1.0)
t0_analytical = 0.5*pi/sqrt(k)

##
f = Figure();

ax = Axis(f[1, 1]; xlabel="t",ylabel="x(t)",title="Harmonic Oscillator")
to_value = x -> x.value
lines!(to_value.(t_array), to_value.(x_array), color=:blue)
f


##
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