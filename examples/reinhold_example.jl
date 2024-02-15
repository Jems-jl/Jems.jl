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
#open("example_options.toml", "w") do file
#    write(file,
#          """
#          [remesh]
#          do_remesh = true
#          delta_log10P_max = 0.1 #1e5 #1e5
#          delta_log10r_max = 0.05 #1e5 #1e-1
#          delta_dm_max = 2e30
#          """)
#end
#StellarModels.set_options!(sm.opt, "./example_options.toml")

open("example_options.toml", "w") do file
    write(file,
          """
          [remesh]
          do_remesh = true
          delta_dm_max = 3.3e30
          delta_log10r_max = 0.015 #1e5 #1e-1
          delta_log10P_max = 0.03 #1e5 #1e5

          [solver]
          newton_max_iter_first_step = 1000
          newton_max_iter = 200

          [timestep]
          dt_max_increase = 10.0
          delta_R_limit = 0.01
          delta_Tc_limit = 0.01

          [termination]
          max_model_number = 50
          max_center_T = 4e7

          [plotting]
          do_plotting = false
          wait_at_termination = false
          plotting_interval = 1

          window_specs = ["HR", "profile", "history"]
          window_layouts = [[1, 1],  # arrangement of plots
                            [2, 1],
                            [3, 1]
                            ]

          profile_xaxis = 'mass'
          profile_yaxes = ['log10_T']
          profile_alt_yaxes = ['X','Y']

          history_xaxis = 'star_age'
          history_yaxes = ['R_surf']
          history_alt_yaxes = ['T_center']

          [io]
          profile_interval = 50
          terminal_header_interval = 100
          terminal_info_interval = 100

          """)
end
StellarModels.set_options!(sm.opt, "./example_options.toml")
rm(sm.opt.io.hdf5_history_filename; force=true)
rm(sm.opt.io.hdf5_profile_filename; force=true)

delta_dm_max = sm.opt.remesh.delta_dm_max
delta_log10r_max = sm.opt.remesh.delta_log10r_max
delta_log10P_max = sm.opt.remesh.delta_log10P_max

print("logP_max -> ")
println(delta_log10P_max)
print("logr_max -> ")
println(delta_log10r_max)
print("dm_max -> ")
println(delta_dm_max)
##

using GLMakie


println("nCells = ", sm.nz)        

f = Figure(size=(600, 600))
idx = 1:sm.nz-1
idx2 = 1:sm.nz
dm = [sm.dm[i] for i in 1:sm.nz-1]
lnr = [sm.psi.lnr[i] for i in 1:sm.nz]/log(10)
lnP = [sm.psi.lnP[i] for i in 1:sm.nz]/log(10)
dlnr = lnr[2:sm.nz] - lnr[1:sm.nz-1]
dlnP = - (lnP[2:sm.nz] - lnP[1:sm.nz-1])
ax1 = Axis(f[1,1], xlabel=L"\mathrm{idx}", ylabel=L"dM", yticklabelcolor=:blue)
ax21 = Axis(f[2,1], xlabel=L"\mathrm{idx}", yticklabelcolor=:blue, ylabel=L"\Delta \log_{10}r")
ax22 = Axis(f[2,1], xlabel=L"\mathrm{idx}", yticklabelcolor=:red, yaxisposition=:right, ylabel=L"\log_{10}r")
ax31 = Axis(f[3,1], xlabel=L"\mathrm{idx}", yticklabelcolor=:blue, ylabel=L"\log_{10}P")
ax32 = Axis(f[3,1], xlabel=L"\mathrm{idx}", yticklabelcolor=:red, yaxisposition=:right, ylabel=L"\log_{10}P")
hidespines!(ax22)
hidexdecorations!(ax22)
hidespines!(ax32)
hidexdecorations!(ax32)
lines!(ax1, idx, dm, color=:blue)
lines!(ax1, idx, delta_dm_max*ones(idx), linestyle=:dash, color=:blue)
lines!(ax21, idx, dlnr, color=:blue)
lines!(ax21, idx,  delta_log10r_max*ones(idx), linestyle=:dash, color=:blue)
lines!(ax22, idx2, lnr, color=:red) 
lines!(ax31, idx, dlnP, color=:blue)
lines!(ax31, idx, delta_log10P_max*ones(idx), linestyle=:dash, color=:blue)
lines!(ax32, idx2, lnP, color=:red) 
f

##

rm(sm.opt.io.hdf5_history_filename; force=true)
rm(sm.opt.io.hdf5_profile_filename; force=true)
sm = StellarModel(varnames, structure_equations, nz, nextra,
                  remesh_split_functions, net, eos, opacity);
StellarModels.set_options!(sm.opt, "./example_options.toml")
StellarModels.n_polytrope_initial_condition!(n, sm, MSUN, 100 * RSUN; initial_dt=10 * SECYEAR)
#StellarModels.remesher!(sm)
sm = Evolution.do_evolution_loop!(sm);

##
##
 



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
