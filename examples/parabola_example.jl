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
