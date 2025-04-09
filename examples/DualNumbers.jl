#=
# DualSupport.jl

An important part of the design of Jems is the thorough use of Dual numbers to compute partial derivatives with respect to all independent variables across the entire code. By independent variables traditionally in a Lagrangian stellar evolution code one refers to two thermodynamic properties (eg $T$ and $\rho$), the luminosity $L$ and the radius $r$ which are functions of the mass coordinate, plus the mass fractions of all isotopes that are being tracked. This means we have a total number of independent variables $n_\mathrm{vars}=4+n_\mathrm{species}$, where $n_\mathrm{species}$ is the number of isotopes included in the model. It can be hard to track how these utilities are used within the code, so the idea of this example is to show how they can be used for a much simpler physical system.

### Dual numbers and caches

Let's first illustrate some basics about Dual numbers. We can create a Dual number using the implementation of the package `ForwardDiff.jl`. We consider functions of two variables $x$ and $y$ for which we will want partial derivatives. We can initialize dual numbers for these as (TODO: maybe show usage of tags to avoid autodiff confusion)
=#

import ForwardDiff: Dual, Tag

xfloat = 3.0
yfloat = 2.0
x = Dual(xfloat,(1.0,0.0))
y = Dual(yfloat,(0.0,1.0))

##
#=
In here the first number given to the `Dual` constructor is the value of the variable, meaning I set $x=3$ and $y=3$, while the tuple that follows are the partial derivatives with respect to $x$ and $y$. For each of these we just have $\partial y/\partial y=1$ and $\partial x/\partial x=1$. Next we can perform some operations with these numbers and verify that we obtain the expected result and partial derivatives (assert will throw an exception if the approximate equality is not satisfied):
=#

function test(x,y)
    #the partial derivatives for this one are
    #∂f/∂x = cos(y), ∂f/∂y=-x*sin(y)
    return x*cos(y)
end
result = test(x,y)

@assert result.value ≈ test(xfloat,yfloat)
@assert result.partials[1] ≈ cos(yfloat)
@assert result.partials[2] ≈ -xfloat*sin(yfloat)

##
#=
One important issue with the use of dual numbers is that constantly creating them can produce a lot of allocations leading to significant slowdown. We can work around this issue by using a cache, which is precisely what the struct `StarDiffCache` is meant to do (the implementation of this struct was adapted from `PreallocationTools.jl`). Considering a case with $n_\mathrm{vars}=2$, we need a cache that can store three values (that of the variable and that of the partial derivatives). We make here a cache for both $x$ and $y$, this is very low level and very likely you will not need to interact with it.
=#

import Jems.DualSupport: StarDiffCache, get_dual

cache_x = StarDiffCache{3,Float64}(zeros(3))
cache_x.dual_data[2] = 1.0 # dual_data[1] contains the value, dual_data[2] is ∂x/∂x=1, dual_data[3] is ∂x/∂y=0

cache_y = StarDiffCache{3,Float64}(zeros(3))
cache_y.dual_data[3] = 1.0; # dual_data[1] contains the value, dual_data[2] is ∂x/∂x=1, dual_data[3] is ∂x/∂y=0

#we can repeat the example from above using these
#first set the value of each dual number, which goes into the first entry of `dual_data`
cache_x.dual_data[1] = xfloat
cache_y.dual_data[1] = yfloat

x = get_dual(cache_x)
y = get_dual(cache_y)
result = test(x,y)

@assert result.value ≈ test(xfloat,yfloat)
@assert result.partials[1] ≈ cos(yfloat)
@assert result.partials[2] ≈ -xfloat*sin(yfloat)

##
#=
### `CellDualData` for calculations in a three-point stencil

An extra level of complexity is concerned with how the equations of stellar structure and evolution are solved. Consider for instance the equation of hydrostatic equilibrium:

$$\frac{\partial P}{\partial m}=-\frac{Gm}{4\pi r^4}$$

Although we are looking for continuous solutions $P(m)$ and $r(m)$, for numerical purposes we need a discretization. The precise locations where pressure and radii are defined for a spherical shell can differ (a staggered mesh, see for instance Figure 9 of Paxton et al. 2011). For simplicity lets define all properties at cell centers, consider all cells of equal mass $\Delta m$ and estimate derivatives with a three point stencil. The above equation can then be written as:

$$\frac{P_{i+1}-P_{i-1}}{2\Delta m} + \frac{G m_i}{4\pi r_i^4} = 0$$

Here $P_{i}$, $\r_{i}$ and $m_{i}$ correspond to the pressure, radius and mass at the center of cell $i$ (where the index increases towards the center). In a Lagrangian code we have fixed cells, meaning $m_i$ is known, and we want to find the values for all other variables in all zones (in this case $P_i$ and $r_i$) for which the above equation is fulfilled. For clarity in this example, lets set everything except pressure and radii to unity,

$$P_{i+1}-P_{i-1}+r_i^{-4}=0.$$

The idea is that we will start with a guess for the solution to this equation (in an actual evolutionary model this can be the previous step properties or an initial condition defined with, for example, a solution to the Lane-Emden equation). The guess contains all values for pressure and density everywhere in the model, but since they need not satisfy hydrostatic equilibrium exactly we can define a residual to the equation of hydrostatic equilibrium at this cell,

$$f_i = P_{i+1}-P_{i-1}+r_i^{-4}.$$

In reality we have more differential equations of stellar structure that we need to consider (normally one per independent variable), and to solve them we make use of a Newton-Rhapson solver which attempts to find a set of values for the independent variables across the entire model that satisfy all the discretized equations. For this purpose we need not only evaluate the residual, but also evaluate all relevant non-zero partial derivatives, which in this case are

$$\frac{\partial f_i}{\partial P_{i+1}}=1,\quad \frac{\partial f_i}{\partial P_{i-1}}=-1,\quad \frac{\partial f_i}{\partial r_i}=-4 r_i^{-5}.$$

So from this we see, we are considering two independent variables $P$ and $r$, but the partial derivatives we need from the residuals need to be taken against the independent variables above, below and at the present cell. In general, for each residual we need $3n_\mathrm{vars}$ partial derivatives, although typically many are equal to zero. A standard practice is just to hard code the partial derivatives as determined from the analytical expressions, but for complex calculations this can be very cumbersome and error prone. So the idea is to setup automatic differentiation tools that take care of this. This is the purpose of the `CellDualData` struct. We start by initializing 6 different `CellDualData` to represent densities and pressures around a point in the three-point stencil:
=#

using Jems.DualSupport

nvars = 2

P_m1 = CellDualData(nvars, Float64; is_ind_var=true, ind_var_i=1) # this is pressure at the cell below
P_00 = CellDualData(nvars, Float64; is_ind_var=true, ind_var_i=1) # this is pressure at the current cell
P_p1 = CellDualData(nvars, Float64; is_ind_var=true, ind_var_i=1) # this is pressure at the cell above
#initialize them with arbitrary values
update_cell_dual_data_value!(P_m1, 0.9)
update_cell_dual_data_value!(P_00, 1.0)
update_cell_dual_data_value!(P_p1, 1.0)

#same for density
r_m1 = CellDualData(nvars, Float64; is_ind_var=true, ind_var_i=2) # this is pressure at the cell below
r_00 = CellDualData(nvars, Float64; is_ind_var=true, ind_var_i=2) # this is pressure at the current cell
r_p1 = CellDualData(nvars, Float64; is_ind_var=true, ind_var_i=2) # this is pressure at the cell above
#initialize them with arbitrary values
update_cell_dual_data_value!(r_m1, 1.0)
update_cell_dual_data_value!(r_00, 1.0)
update_cell_dual_data_value!(r_p1, 1.0);

##
#=
With this in place we can create dual numbers and operate on them. We do this using the `get_00_dual`, `get_m1_dual` and `get_p1_dual`. This ensures the partial derivatives are filled in the way we want. To illustrate this, let's see the partials produced by the following
=#

P_m1_dual = get_m1_dual(P_m1)
r_m1_dual = get_m1_dual(r_m1)
P_00_dual = get_00_dual(P_00)
r_00_dual = get_00_dual(r_00)
P_p1_dual = get_p1_dual(P_p1)
r_p1_dual = get_p1_dual(r_p1)

@show P_m1_dual.partials
@show r_m1_dual.partials
@show P_00_dual.partials
@show r_00_dual.partials
@show P_p1_dual.partials
@show r_p1_dual.partials;

##
#=
As you can see we have dual numbers with a total of 6 partial derivatives. When combining these together to compute a residual $f_i$, the partial derivatives will contain (in order)

$$\frac{\partial f_i}{\partial P_{i-1}},\;\frac{\partial f_i}{\partial r_{i-1}},\;\frac{\partial f_i}{\partial P_{i}},\;\frac{\partial f_i}{\partial r_{i}},\;\frac{\partial f_i}{\partial P_{i+1}},\;\frac{\partial f_i}{\partial r_{i+1}}.$$

We can directly evaluate the equation for $f_i$ that we used above and verify this works.
=#

f_i = P_p1_dual-P_m1_dual+r_00_dual^(-4)

#verify non-zero partials
@assert f_i.partials[1] ≈ -1 # This is testing \partial f_i/ \partial P_{i-1}
@assert f_i.partials[4] ≈ -4*r_00_dual.value^5 # This is testing \partial f_i/ \partial r_{i}
@assert f_i.partials[5] ≈ 1 # This is testing \partial f_i/ \partial P_{i}

##
#=

### Evaluating and storing results into `CellDualData` instances (TODO)

### A simple example including a Newton_Rhapson solver (TODO)

Here I want to show a simple system with a known solution, a vertical tube with an ideal gas fixed at constant temperatur, where the idea is to determine the density profile at arbitrary heights when there is a constant gravity.

### Further complications, `FaceDualData` (TODO)

Properties evaluated at faces are defined differently, need to also summarize why and how.

=#