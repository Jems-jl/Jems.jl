import ForwardDiff.Dual

#defining differential equations
dy_dt(y,vy)  = vy
dvy_dt(y,vy) = -9.81 #m/s²

max_nb_steps = 100_000

# defining initial conditions as dual numbers
t0      = Dual(0.0  ,0.0)
Delta_t = Dual(0.001,0.0)
y0 = h  = Dual(10.0 ,1.0) # we keep track of derivatives with respect to height
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