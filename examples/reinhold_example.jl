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
StellarModels.update_stellar_model_properties!(sm)
Evolution.eval_jacobian_eqs!(sm)

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
          delta_log10P_max = 100.0 #1eh
          #delta_log10r_max = 1e-1 #1e9
          """)
end
StellarModels.set_options!(sm.opt, "./example_options.toml")



##

using CairoMakie
using LaTeXStrings
using Printf
 
# using printf macro with @ sign 
StellarModels.remesher!(sm)
println(sm.nz)        
##
#[abs(sm.psi.lnr[i] - sm.psi.lnr[i+1])/log(10) for i in 1:sm.nz]
#maxval, idx = findmax()
println( maxval) 
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