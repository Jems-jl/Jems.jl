using BenchmarkTools
using Jems.Chem
using Jems.Constants
using Jems.EOS
using Jems.Opacity
using Jems.NuclearNetworks
using Jems.Turbulence
using Jems.StellarModels
using Jems.Evolution
using Jems.ReactionRates

##

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4], [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                    nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);


n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
StellarModels.cycle_props!(sm);
StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

a = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end

BenchmarkTools.save("BM1_Toy_2ISO_A.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Toy_2ISO_B.json", b)


##

# 3 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :D2, :He4], [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                    nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);


n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
StellarModels.cycle_props!(sm);
StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

a = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end

BenchmarkTools.save("BM1_Toy_3ISO_A.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Toy_3ISO_B.json", b)

##

# 4 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :D2, :He3, :He4], [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                    nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);


n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
StellarModels.cycle_props!(sm);
StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

a = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end

BenchmarkTools.save("BM1_Toy_4ISO_A.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Toy_4ISO_B.json", b)

##


# 5 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :D2, :He3, :Li3, :He4], [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                    nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);


n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
StellarModels.cycle_props!(sm);
StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

a = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end

BenchmarkTools.save("BM1_Toy_5ISO_A.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Toy_5ISO_B.json", b)

##

# 6 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :D2, :T3, :He3, :Li3,
                    :He4],
                    [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                    nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);


n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
StellarModels.cycle_props!(sm);
StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

a = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end

BenchmarkTools.save("BM1_Toy_6ISO_A.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Toy_6ISO_B.json", b)

##


# 7 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :D2, :T3, :He3, :Li3,
                    :He4, :Li7],
                    [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                    nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);


n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
StellarModels.cycle_props!(sm);
StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

a = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end

BenchmarkTools.save("BM1_Toy_7ISO_A.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Toy_7ISO_B.json", b)

##

# 8 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :D2, :T3, :He3, :Li3,
                    :He4, :Li7, :Be9],
                    [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                    nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);


n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
StellarModels.cycle_props!(sm);
StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

a = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end

BenchmarkTools.save("BM1_Toy_8ISO_A.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Toy_8ISO_B.json", b)

##


# 9 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :D2, :T3, :He3, :Li3,
                    :He4, :Li7, :Be9, :B11],
                    [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                    nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);


n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
StellarModels.cycle_props!(sm);
StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

a = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end

BenchmarkTools.save("BM1_Toy_9ISO_A.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Toy_9ISO_B.json", b)

##

# 10 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :D2, :T3, :He3, :Li3,
                    :He4, :Li7, :Be9, :B11, :C12],
                    [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                    nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);


n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
StellarModels.cycle_props!(sm);
StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

a = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end

BenchmarkTools.save("BM1_Toy_10ISO_A.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Toy_10ISO_B.json", b)

##

# 15 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :D2, :T3, :He3, :Li3,
                    :He4, :Li7, :Be9, :B11, :C12,
                    :N14, :O16, :F19, :Ne20, :Na23,],
                    [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                    nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);


n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
StellarModels.cycle_props!(sm);
StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

a = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end

BenchmarkTools.save("BM1_Toy_15ISO_A.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Toy_15ISO_B.json", b)

##

# 20 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :D2, :T3, :He3, :Li3,
                    :He4, :Li7, :Be9, :B11, :C12,
                    :N14, :O16, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32],
                    [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                    nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);


n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
StellarModels.cycle_props!(sm);
StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

a = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end

BenchmarkTools.save("BM1_Toy_20ISO_A.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Toy_20ISO_B.json", b)

##

# 25 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :D2, :T3, :He3, :Li3,
                    :He4, :Li7, :Be9, :B11, :C12,
                    :N14, :O16, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45],
                    [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                    nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);


n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
StellarModels.cycle_props!(sm);
StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

a = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end

BenchmarkTools.save("BM1_Toy_25ISO_A.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Toy_25ISO_B.json", b)

##

# 30 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :D2, :T3, :He3, :Li3,
                    :He4, :Li7, :Be9, :B11, :C12,
                    :N14, :O16, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    :Ti48, :V51, :Cr52, :Mn55, :Fe56],
                    [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                    nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);


n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
StellarModels.cycle_props!(sm);
StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

a = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end

BenchmarkTools.save("BM1_Toy_30ISO_A.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Toy_30ISO_B.json", b)

##


# 35 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :D2, :T3, :He3, :Li3,
                    :He4, :Li7, :Be9, :B11, :C12,
                    :N14, :O16, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    :Co59, :Ni58, :Cu63, :Zn64, :Ga69],
                    [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                    nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);


n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
StellarModels.cycle_props!(sm);
StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

a = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end

BenchmarkTools.save("BM1_Toy_35ISO_A.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Toy_35ISO_B.json", b)


##


# 40 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :D2, :T3, :He3, :Li3,
                    :He4, :Li7, :Be9, :B11, :C12,
                    :N14, :O16, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    :Ge74, :As75, :Se80, :Br79, :Kr84],
                    [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                    nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);


n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
StellarModels.cycle_props!(sm);
StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

a = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end

BenchmarkTools.save("BM1_Toy_40ISO_A.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Toy_40ISO_B.json", b)

##


# 45 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :D2, :T3, :He3, :Li3,
                    :He4, :Li7, :Be9, :B11, :C12,
                    :N14, :O16, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    :Ge74, :As75, :Se80, :Br79, :Kr84,
                    :Rb85, :Sr88, :Y89, :Zr90, :Nb93],
                    [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                    nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);


n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
StellarModels.cycle_props!(sm);
StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

a = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end

BenchmarkTools.save("BM1_Toy_45ISO_A.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Toy_45ISO_B.json", b)

##

# 50 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :D2, :T3, :He3, :Li3,
                    :He4, :Li7, :Be9, :B11, :C12,
                    :N14, :O16, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    :Ge74, :As75, :Se80, :Br79, :Kr84,
                    :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:toy_rates, :toy_pp), (:toy_rates, :toy_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(true)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                    nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);


n = 3
StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0, Chem.abundance_lists[:ASG_09], MSUN,
                                             100 * RSUN; initial_dt=10 * SECYEAR)
StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
StellarModels.cycle_props!(sm);
StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

a = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end

BenchmarkTools.save("BM1_Toy_50ISO_A.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Toy_50ISO_B.json", b)


##





