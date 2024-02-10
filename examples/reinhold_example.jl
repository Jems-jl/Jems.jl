#=
# NuclearBurning.jl

This notebook provides a simple example of a star with simplified microphysics undergoing nuclear burning.
Import all necessary Jems modules. We will also do some benchmarks, so we import BenchmarkTools as well.
=#
using BenchmarkTools
using Jems.Chem
using Jems.Constants
using Jems.EOS
using Jems.Opacity
using Jems.NuclearNetworks
using Jems.StellarModels
using Jems.Evolution
using Jems.ReactionRates

using CairoMakie
#using LaTeXStrings
using Printf
 
##
#=
### Model creation

We start by creating the stellar model. In this example we consider a model with 6 independent variables, two of which
correspond to composition. The independent variables here are $\ln(P)$, $\ln(T)$, $\ln(r)$, the luminosity $L$ and the
mass fractions of Hydrogen and Helium.

The Evolution module has pre-defined equations corresponding to these variables, which we provide here. For now, only a
simple (fully ionized) ideal gas law EOS is available. Similarly, only a simple simple electron scattering opacity equal
to $\kappa=0.2(1+X)\;[\mathrm{cm^2\;g^{-1}}]$ is available.
=#

varnames = [:lnρ, :lnT, :lnr, :lum]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1,:He4,:C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(false)
opacity = Opacity.SimpleElectronScatteringOpacity()
sm = StellarModel(varnames, structure_equations, nz, nextra,
                  remesh_split_functions, net, eos, opacity);

n=3
StellarModels.n_polytrope_initial_condition!(n, sm, MSUN, 100 * RSUN; initial_dt=10 * SECYEAR)
Evolution.set_step_info!(sm, sm.esi)
Evolution.cycle_step_info!(sm);
Evolution.set_step_info!(sm, sm.ssi)
StellarModels.update_stellar_model_properties!(sm, sm.props)
Evolution.eval_jacobian_eqs!(sm)

#println(sm.psi.lnP)

for i in 1:sm.psi.nz
    println(sm.psi.lnP[i])
end
sm.nz

##
#=
### Evolving our model

We can now evolve our star! We will initiate a $1M_\odot$ star with a radius of $100R_\odot$ using an n=1 polytrope (it
would be much better to use n=3 or n=3/2 polytropes, for now I only use this because there is a simple analytical
solution). The star is expected to contract until it ignites hydrogen. We set a few options for the simulation with a
toml file, which we generate dynamically. These simulation should complete in about a thousand steps once it reaches the
`max_center_T` limit.

Output is stored in HDF5 files, and easy to use functions are provided with the Evolution module to turn these HDF5
files into DataFrame objects. HDF5 output is compressed by default.
=#
open("example_options.toml", "w") do file
    write(file,
          """
          [remesh]
          do_remesh = true
          delta_log10P_max = 1e5 #1e5
          delta_log10r_max = 1e5 #1e-1
          delta_dm_max = 1e29
          """)
end
StellarModels.set_options!(sm.opt, "./example_options.toml")

##


delta_log10P_max = sm.opt.remesh.delta_log10P_max
delta_log10r_max = sm.opt.remesh.delta_log10r_max
delta_dm_max = sm.opt.remesh.delta_dm_max

print("logP_max -> ")
println(delta_log10P_max)
print("logr_max -> ")
println(delta_log10r_max)
print("dm_max -> ")
println(delta_dm_max)
##



# we first do cell splitting
#do_split = Vector{Bool}(undef, sm.nz) 
#do_split .= false
###maxPinstar = 0
###maxrinstar = 0
#for i in 1:sm.nz-1
#    delta_log10P = abs(sm.psi.lnP[i] - sm.psi.lnP[i+1])/log(10)
#    #if delta_log10P > maxPinstar
#    #    maxPinstar = delta_log10P
#    delta_log10P_max = sm.opt.remesh.delta_log10P_max
#    delta_log10r = abs(sm.psi.lnr[i] - sm.psi.lnr[i+1])/log(10)
#    #if delta_log10r > maxrinstar
#    #    maxrinstar = delta_log10r
#    delta_log10r_max = sm.opt.remesh.delta_log10r_max
#    if (delta_log10r > delta_log10r_max)
#        println(i, ", logr")
#        # if the condition is satisfied, we split the largest of the two cells
#        if sm.dm[i] > sm.dm[i+1]
#            do_split[i] = true
#        else
#            do_split[i+1] = true
#        end
#    elseif (delta_log10P > delta_log10P_max) 
#        println(i, ", logP")
#        # if the condition is satisfied, we split the largest of the two cells
#        if sm.dm[i] > sm.dm[i+1]
#            do_split[i] = true
#        else
#            do_split[i+1] = true
#        end
#    end
#end
### we are ignoring the edges for now
### do_split[sm.nz] = false
#extra_cells = sum(do_split)
# if allocated space is not enough, we need to reallocate everything
#if sm.nz + extra_cells > length(sm.dm)
#    sm = adjusted_stellar_model_data(sm, sm.nz + extra_cells, sm.nextra);
#end
#for i=sm.nz:-1:2
#    if extra_cells == 0
#        break
#    end
#    if !do_split[i]
#        # move everything upwards by extra_cells
#        for j in 1:sm.nvars
#            sm.ind_vars[(i+extra_cells-1)*sm.nvars+j] = sm.ind_vars[(i-1)*sm.nvars+j]
#        end
#        # RTW: why do we need these lines here? Is mass not part of the ind_vars?
#        sm.m[i+extra_cells] = sm.m[i]
#        sm.dm[i+extra_cells] = sm.dm[i]
#        if i > sm.nz - 10
#            println("i, extra = ", i, " ", extra_cells)
#            println(sm.m[i])
#            println(sm.m[i+extra_cells])
#            println(sm.dm[i])
#            println(sm.dm[i+extra_cells])
#            println()
#        end
#    else
#        println("This index is getting split ", i)
#        println("Extra cells = ", extra_cells)
#        #split the cell and lower extra_cells by 1
#        #dm_00 = sm.dm[i]
#        #var_00 = view(sm.ind_vars, ((i-1)*sm.nvars + 1):((i-1)*sm.nvars + sm.nvars))
#        #if i > 1
#        #    dm_m1 = sm.dm[i-1]
#        #    var_m1 = view(sm.ind_vars, ((i-2)*sm.nvars + 1):((i-2)*sm.nvars + sm.nvars))
#        #else
#        #    dm_m1 = NaN
#        #    var_m1 = []
#        #end
#        #if i < sm.nz
#        #    dm_p1 = sm.dm[i+1]
#        #    var_p1 = view(sm.ind_vars, ((i)*sm.nvars + 1):((i)*sm.nvars + sm.nvars))
#        #else
#        #    dm_p1 = NaN
#        #    var_p1 = []
#        #end
#        #varnew_low = view(sm.ind_vars, 
#        #                  ((i+extra_cells-2)*sm.nvars+1):((i+extra_cells-2)*sm.nvars+sm.nvars))
#        #varnew_up = view(sm.ind_vars, 
#        #                  ((i+extra_cells-1)*sm.nvars+1):((i+extra_cells-1)*sm.nvars+sm.nvars))
#
#        #for remesh_split_function in sm.remesh_split_functions
#        #    remesh_split_function(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1,
#        #                            varnew_low, varnew_up)
#        #end
#
#        sm.m[i+extra_cells] = 0.5*sm.m[i]
#        sm.m[i+extra_cells-1] = 0.5*sm.m[i]
#        sm.dm[i+extra_cells] = 0.5*sm.dm[i]
#        sm.dm[i+extra_cells-1] = 0.5*sm.dm[i]
#
#        extra_cells = extra_cells - 1
#    end
#end
#sm.nz = sm.nz + sum(do_split)


# using printf macro with @ sign 
#println("nCells = ", sm.nz)        

StellarModels.remesher!(sm)

#print("Sum of split is: ")
#println(sum(do_split))
println("nCells = ", sm.nz)        

f = Figure()
idx = 1:sm.nz-1
#delta_log10_P = [abs(sm.psi.lnP[i] - sm.psi.lnP[i+1])/log(10) for i in 1:sm.nz-1]
delta_dm = [abs(sm.psi.dm[i] - sm.psi.dm[i+1]) for i in 1:sm.nz-1]
ax = Axis(f[1,1], xlabel=L"\mathrm{idx}", ylabel=L"\mathrm{dM}")
lines!(ax, idx, delta_dm)
f






#StellarModels.profile_output_options["log10_P"][2](sm, k) for k in 1:sm.nz]



#for i in 1:sm.nz-1
#    #delta_log10P = abs(sm.psi.lnP[i] - sm.psi.lnP[i+1])/log(10)
#    delta_log10r = abs(sm.psi.lnr[i] - sm.psi.lnr[i+1])/log(10)
#    
#    #if (delta_log10P > delta_log10P_max) 
#    if (delta_log10r > delta_log10r_max) 
#        print(i)
#        print(" => ")
#        #println(delta_log10P)
#        println(delta_log10r)
#    end
#end
##


# we first do cell splitting
#do_split = Vector{Bool}(undef, sm.nz) 
#do_split .= false
#maxPinstar = 0
#maxrinstar = 0

##
    #if delta_log10P > maxPinstar
    #    maxPinstar = delta_log10P
    #if delta_log10r > maxrinstar
    #    maxrinstar = delta_log10r
    delta_log10r_max = sm.opt.remesh.delta_log10r_max
    if (delta_log10r > delta_log10r_max)
        # if the condition is satisfied, we split the largest of the two cells
        if sm.dm[i] > sm.dm[i+1]
            do_split[i] = true
        else
            do_split[i+1] = true
        end
    elseif (delta_log10P > delta_log10P_max) 
        # if the condition is satisfied, we split the largest of the two cells
        if sm.dm[i] > sm.dm[i+1]
            do_split[i] = true
        else
            do_split[i+1] = true
        end
    end
end












##
#[abs(sm.psi.lnr[i] - sm.psi.lnr[i+1])/log(10) for i in 1:sm.nz]
#maxval, idx = findmax()
#println( maxval) 
maxval, idx = findmax([abs(sm.psi.lnP[i] - sm.psi.lnP[i+1])/log(10) for i in 1:sm.nz])
println( maxval) 


##

idx = 1:sm.nz-1
println(last(idx))
println(idx[1])
log10_P = [StellarModels.profile_output_options["log10_P"][2](sm, k) for k in 1:sm.nz]
log10_r = [StellarModels.profile_output_options["log10_r"][2](sm, k) for k in 1:sm.nz]
delta_log10_P = log10_P[2:end] - log10_P[1:end-1]
delta_log10_r = log10_r[2:end] - log10_r[1:end-1]
#Psize = size(delta_log10_P)
#println(Psize)

#mass = [StellarModels.profile_output_options["mass"][2](sm, k) for k in 1:sm.nz]
#deltaP = log10_P[k+1]-log10_P[k]

f = Figure()
ax = Axis(f[1,1], xlabel=L"\mathrm{idx}", ylabel=L"\mathrm{P}\;[P_\odot?]")
lines!(ax, idx, delta_log10_P)
lines!(ax, idx, delta_log10_r)
f


##
iρ = sm.vari[:lnρ]
lnρ = [sm.ind_vars[ (k-1)*sm.nvars+iρ] for k in 1:sm.nz]
ir = sm.vari[:lnr]
lnr = [sm.ind_vars[(k-1)*sm.nvars+ir] for k in 1:sm.nz]

f = Figure()
ax = Axis(f[1,1], xlabel=L"\mathrm{mass}\;[M_\odot]", ylabel=L"\mathrm{T}\;[T_\odot]")
lines!(ax, lnr, lnρ)
f









##
#=
### Initialize StellarModel and evaluate equations and jacobian

We do not have a working initial condition yet. We require pressure, temperature profiles. One simple available initial
condition is that of an n=0 polytrope. This sets the pressure and density and computes the temperature from the EOS. The
luminosity is initialized by assuming pure radiative transport for the temperature gradient produced by the polytrope.

The normal evolution loop will store the information at the end of the step into an attribute of type `StellarStepInfo`,
stored at `sm.esi` (_end step info_). After initializing our polytrope we can mimic that behavior by calling 
`set_end_step_info!(sm)`. We then 'cycle' this info into the information of a hypothetical previous step with
`cycle_step_info`, so now `sm.psi` contains our initial condition. Finally we call `set_start_step_info` to use `sm.psi`
(_previous step info_) to populate the information needed before the Newton solver in `sm.ssi` (_start step info_).
At last we are in position to evaluate the equations and compute the Jacobian.
=#

##
