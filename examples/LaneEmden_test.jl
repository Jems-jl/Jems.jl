5<10<20<=300
##
println("Hello World!")

xvalues = [0,10,20]
yvalues = [0,5,10]

function linear_interpolation(xvalues, yvalues)
    function θ_n(x)
        for i in eachindex(xvalues)
            if xvalues[i] <= x < xvalues[i+1]
                return yvalues[i] + (yvalues[i+1] - yvalues[i]) / (xvalues[i+1] - xvalues[i]) * (x - xvalues[i])
            end
        end
    end
    return θ_n
end
linear_interpolation(xvalues, yvalues)(1)
myfunction = linear_interpolation(xvalues, yvalues)
myfunction(18)


##
println("Start Run")
dydx(x,y,z,n) = z #derivative dy/dx = z
dzdx(x,y,z,n) = -y^n -2*z/x #derivative dz/dx = -y^n -2*z/x
function RK_step(x,y,z,n,Δx)
    #this function returns the new (y,z) values, given the old values(x,y,z)
    #preparing all constants for the Runge-Kutta method
    k₁ = Δx*dydx(x,y,z,n)
    l₁ = Δx*dzdx(x,y,z,n)
    k₂ = Δx*dydx(x+Δx/2,y+k₁/2,z+l₁/2,n)
    l₂ = Δx*dzdx(x+Δx/2,y+k₁/2,z+l₁/2,n)
    k₃ = Δx*dydx(x+Δx/2,y+k₂/2,z+l₂/2,n)
    l₃ = Δx*dzdx(x+Δx/2,y+k₂/2,z+l₂/2,n)
    k₄ = Δx*dydx(x+Δx,y+k₃,z+l₃,n)
    l₄ = Δx*dzdx(x+Δx,y+k₃,z+l₃,n)
    return (
        y+k₁/6+k₂/3+k₃/3+k₄/6, #new y value
        z+l₁/6+l₂/3+l₃/3+l₄/6, #new z value
        )
end
#approximation of y and z for small ξ
y_smallx(x,n) = 1 - 1/6*x^2 + n/120*x^4 -n*(8*n-5)/1520*x^6
z_smallx(x,n) = - 1/3*x + n/30*x^3 -3*n*(8*n-5)/760*x^5;


#set up grid in x 
Δx = 1e-4
n = 0
println("n=$n")
nsteps = 200000 #putting a maximum on the number of steps
xvals = LinRange(Δx,nsteps*Δx,nsteps)

#initialize first value of y and z using series approximation (i=1)
yvals = zeros(nsteps); zvals = zeros(nsteps)
yvals[1] = y_smallx(Δx,n); zvals[1] = z_smallx(Δx,n)
#print type of yvals

 

function endOfLoop!(xrange::LinRange, yrange::Vector{Float64}, zrange::Vector{Float64},endIndex::Int64)
    println("END OF LOOP ACTIVATED, i = $endIndex")
    slope = (yrange[endIndex-1] - yrange[endIndex-2]) / (xrange[endIndex-1] - xrange[endIndex-2])
    xlastoriginal = xrange[endIndex]
    xlast = xrange[endIndex-1] - yrange[endIndex-1] / slope
    newxvals = zeros(endIndex)
    newxvals[1:endIndex-1] = xrange[1:endIndex-1]
    newxvals[endIndex] = xlast #linearly interpolated xvalue, this last value is $ξ_1$
    yvals[endIndex] = 0.0 #manually set this value to zero
    newyvals = yrange[1:endIndex]; newzvals = zrange[1:endIndex] #we truncate the arrays
    #return (newxvals, newyvals, newzvals)
    return (newxvals, push!(yrange[1:endIndex-1],0.0))
end

bla = xvals[1]



#perform Runge-Kutta integration
for i in 2:nsteps
    x = xvals[i-1]; 
    y = yvals[i-1]; z = zvals[i-1]
    #calculating the RK constants
    k₁ = Δx*dydx(x,y,z,n)
    l₁ = Δx*dzdx(x,y,z,n)
    ynew = y+k₁/2
    if ynew < 0.0
        xvals, yvals = endOfLoop!(xvals, yvals, zvals,i)
        break
    end
    k₂ = Δx*dydx(x+Δx/2,ynew,z+l₁/2,n)
    l₂ = Δx*dzdx(x+Δx/2,ynew,z+l₁/2,n)
    ynew = y+k₂/2
    if ynew < 0.0
        xvals, yvals= endOfLoop!(xvals, yvals, zvals,i)
        break
    end
    k₃ = Δx*dydx(x+Δx/2,ynew,z+l₂/2,n)
    l₃ = Δx*dzdx(x+Δx/2,ynew,z+l₂/2,n)
    ynew = y+k₃
    if ynew < 0.0
        xvals, yvals = endOfLoop!(xvals, yvals, zvals,i)
        break
    end
    k₄ = Δx*dydx(x+Δx,ynew,z+l₃,n)
    l₄ = Δx*dzdx(x+Δx,ynew,z+l₃,n)
    
    yvals[i] = y+k₁/6+k₂/3+k₃/3+k₄/6 #new y value
    zvals[i] = z+l₁/6+l₂/3+l₃/3+l₄/6 #new z value
end  
        


using CairoMakie

println("ξ_1 = $(xvals[end])")
println("y(ξ_1) = $(yvals[end])")
f = Figure()
ax = Axis(f[1, 1])
lines!(ax, xvals, yvals)
scatter!(ax, [xvals[end]], [0.0], markersize = 20) 
println("End Run, show plot")
#x label
ax.xlabel = "ξ"
ax.ylabel = "θ_n(ξ)"
lines!(ax,xvals, 1 .-xvals.^2 ./6)
f
