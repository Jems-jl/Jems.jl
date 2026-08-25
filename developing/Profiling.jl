#=
# NuclearBurning.jl

This notebook provides a simple example of a star with simplified microphysics undergoing nuclear burning.
Import all necessary Jems modules. We will also do some benchmarks, so we import BenchmarkTools as well.
=#
using BenchmarkTools
using PProf
using Profile
using Jems.Chem
using Jems.Constants
using Jems.EOS
using Jems.Opacity
using Jems.NuclearNetworks
using Jems.Turbulence
using Jems.StellarModels
using Jems.Evolution

##
#=
### Model creation
≤
We start by creating the stellar model. In this example we consider a model with 6 independent variables, two of which
correspond to composition. The independent variables here are $\ln(P)$, $\ln(T)$, $\ln(r)$, the luminosity $L$ and the
mass fractions of Hydrogen and Helium.

The Evolution module has pre-defined equations corresponding to these variables, which we provide here. For now, only a
simple (fully ionized) ideal gas law EOS is available. Similarly, only a simple simple electron scattering opacity equal
to $\kappa=0.2(1+X)\;[\text{cm^2\;g^{-1}}]$ is available.
=#

##

net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(StellarModels.DefaultStellarEquationSet(), nz, nextra, net, eos, opacity, turbulence);

##
#=
### Initialize StellarModel and evaluate equations and jacobian

We do not have a working initial condition yet. We require pressure, temperature profiles. One simple available initial
condition is that of an n=1 polytrope. This sets the pressure and density and computes the temperature from the EOS. The
luminosity is initialized by assuming pure radiative transport for the temperature gradient produced by the polytrope.
Information of the model at its present and following step are required at the beginning, the function
`compute_starting_model_properties!` takes care of setting this up.
=#
n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
Evolution.compute_starting_model_properties!(sm)

##
#=
### Benchmarking

The previous code leaves everything ready to solve the linearized system.
For now we make use of a the serial Thomas algorithm for tridiagonal block matrices.
We first show how long it takes to evaluate the Jacobian matrix. This requires two
steps, the first is to evaluate properties across the model (for example, the EOS)
and then evaluate all differential equations and fill the Jacobian. We first benchmark
the evaluation of model properties:
=#
@benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
end

##
#=
And next we benchmark the evaluation of the model equations and construction of the Jacobian:
=#
@benchmark begin
    Evolution.eval_jacobian_eqs!($sm)
end

##
#=

To benchmark the linear solver itself we need to perform
the jacobian evaluation as a setup for the benchmark. This is because the solver
destroys the Jacobian to perform in-place operations.
=#

@benchmark begin
    Evolution.block_tridiagonal_solver!($sm, $sm.solver_data)
end setup=(Evolution.eval_jacobian_eqs!($sm))


## Profiling
PProf.clear()
Profile.Allocs.clear()
Profile.Allocs.@profile begin
    StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
end

PProf.Allocs.pprof()
##

@pprof begin
    Evolution.block_tridiagonal_solver!($sm, $sm.solver_data)
end setup=(Evolution.eval_jacobian_eqs!($sm))
