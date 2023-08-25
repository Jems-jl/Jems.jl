# Evolution

The Evolution module of Jems contains the basic tools needed to combine all other modules to perform stellar evolution.
It allows a fully customizable definition of the equations that are solved together with their boundary conditions.

A model is initialized by specifying its independent variables and equations. The following creates a model with a basic ideal gas EOS and
electron scattering opacity, which is then initialized using an n=1 polytrope of 1$M_\odot$ and $100R_\odot$.
Options for a simulation can be specified with a [toml](https://toml.io/en/) file. Below we directly
create a file with custom options and then load it up into the model. After setting everything up the simulation is run.

```@example
using Jems.Evolution
using Jems.EOS
using Jems.Opacity
using Jems.Constants

nvars = 6
nspecies = 2
varnames = [:lnP, :lnT, :lnr, :lum, :H1, :He4]
structure_equations=[Evolution.equationHSE, Evolution.equationT,
                        Evolution.equationContinuity, Evolution.equationLuminosity,
                        Evolution.equationH1, Evolution.equationHe4]
nz = 1000
eos = EOS.IdealEOS(false)
opacity = Opacity.SimpleElectronScatteringOpacity()
sm = StellarModel(varnames, structure_equations, nvars, nspecies, nz, eos, opacity)

#Initialize the model as n=1 polytrope with an initial timestep of 10 years
Evolution.n1_polytrope_initial_condition(sm, MSUN, 100*RSUN; initial_dt=10*SECYEAR)

#Load custom options
open("options.toml","w") do file
    write(file,"""
            [termination]
            max_model_number = 300
""")
end
Evolution.set_options!(sm.opt, "options.toml")

#run simulation
#Evolution.do_evolution_loop(sm)
```

The results of the simulation are provided in HDF5 format.

## StellarModel.jl

```@autodocs
Modules = [Jems.Evolution]
Pages = ["Evolution/StellarModel.jl"]
```

## Options.jl

```@autodocs
Modules = [Jems.Evolution]
Pages = ["Evolution/Options.jl"]
```

## Equations.jl

```@autodocs
Modules = [Jems.Evolution]
Pages = ["Evolution/Equations.jl"]
```

## Solver.jl

```@autodocs
Modules = [Jems.Evolution]
Pages = ["Evolution/Solver.jl"]
```

## EvolutionLoop.jl

```@autodocs
Modules = [Jems.Evolution]
Pages = ["Evolution/EvolutionLoop.jl"]
```

## InitialCondition.jl

```@autodocs
Modules = [Jems.Evolution]
Pages = ["Evolution/InitialCondition.jl"]
```

## IO.jl

```@autodocs
Modules = [Jems.Evolution]
Pages = ["Evolution/IO.jl"]
```