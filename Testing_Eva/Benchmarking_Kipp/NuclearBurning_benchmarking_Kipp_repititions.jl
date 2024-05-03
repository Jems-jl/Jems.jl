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



# 5 ISOTOPES

println("status check 1")

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16],
                    # :D2, :T3, :He3, :Li3, :Li7,
                    # :Be9, :B11, :F19, :Ne20, :Na23,
                    # :Mg24, :Al27, :Si28, :P31, :S32,
                    # :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    # :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    # :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    # :Ge74, :As75, :Se80, :Br79, :Kr84,
                    # :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    # :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

a1 = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end
BenchmarkTools.save("BM1_Kipp_5ISO_A_V2.json", a1)

# a2 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_5ISO_A_V2_rep2.json", a2)

# a3 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_5ISO_A_V2_rep3.json", a3)


b1 = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end
BenchmarkTools.save("BM1_Kipp_5ISO_B_V2.json", b1)

# b2 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
#     Evolution.thomas_algorithm!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_5ISO_B_V2_rep2.json", b2)

# b3 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
#     Evolution.thomas_algorithm!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_5ISO_B_V2_rep3.json", b3)



# 6 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2], #:T3, :He3, :Li3, :Li7,
                    # :Be9, :B11, :F19, :Ne20, :Na23,
                    # :Mg24, :Al27, :Si28, :P31, :S32,
                    # :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    # :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    # :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    # :Ge74, :As75, :Se80, :Br79, :Kr84,
                    # :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    # :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

a1 = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end
BenchmarkTools.save("BM1_Kipp_6ISO_A_V2.json", a1)

# a2 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_6ISO_A_V2_rep2.json", a2)

# a3 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_6ISO_A_V2_rep3.json", a3)


b1 = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end
BenchmarkTools.save("BM1_Kipp_6ISO_B_V2.json", b1)

# b2 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
#     Evolution.thomas_algorithm!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_6ISO_B_V2_rep2.json", b2)

# b3 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
#     Evolution.thomas_algorithm!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_6ISO_B_V2_rep3.json", b3)


# 7 ISOTOPES
println("status check 2")

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3], # :He3, :Li3, :Li7,
                    # :Be9, :B11, :F19, :Ne20, :Na23,
                    # :Mg24, :Al27, :Si28, :P31, :S32,
                    # :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    # :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    # :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    # :Ge74, :As75, :Se80, :Br79, :Kr84,
                    # :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    # :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

a1 = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end
BenchmarkTools.save("BM1_Kipp_7ISO_A_V2.json", a1)

# a2 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_7ISO_A_V2_rep2.json", a2)

# a3 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_7ISO_A_V2_rep3.json", a3)


b1 = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end
BenchmarkTools.save("BM1_Kipp_7ISO_B_V2.json", b1)

# b2 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
#     Evolution.thomas_algorithm!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_7ISO_B_V2_rep2.json", b2)

# b3 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
#     Evolution.thomas_algorithm!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_7ISO_B_V2_rep3.json", b3)


# 8 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3], # :Li3, :Li7,
                    # :Be9, :B11, :F19, :Ne20, :Na23,
                    # :Mg24, :Al27, :Si28, :P31, :S32,
                    # :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    # :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    # :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    # :Ge74, :As75, :Se80, :Br79, :Kr84,
                    # :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    # :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

a1 = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end
BenchmarkTools.save("BM1_Kipp_8ISO_A_V2.json", a1)

# a2 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_8ISO_A_V2_rep2.json", a2)

# a3 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_8ISO_A_V2_rep3.json", a3)


b1 = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end
BenchmarkTools.save("BM1_Kipp_8ISO_B_V2.json", b1)

# b2 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
#     Evolution.thomas_algorithm!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_8ISO_B_V2_rep2.json", b2)

# b3 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
#     Evolution.thomas_algorithm!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_8ISO_B_V2_rep3.json", b3)

# 9 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3], # :Li7,
                    # :Be9, :B11, :F19, :Ne20, :Na23,
                    # :Mg24, :Al27, :Si28, :P31, :S32,
                    # :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    # :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    # :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    # :Ge74, :As75, :Se80, :Br79, :Kr84,
                    # :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    # :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

a1 = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end
BenchmarkTools.save("BM1_Kipp_9ISO_A_V2.json", a1)

# a2 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_9ISO_A_V2_rep2.json", a2)

# a3 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_9ISO_A_V2_rep3.json", a3)


b1 = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end
BenchmarkTools.save("BM1_Kipp_9ISO_B_V2.json", b1)

# b2 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
#     Evolution.thomas_algorithm!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_9ISO_B_V2_rep2.json", b2)

# b3 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
#     Evolution.thomas_algorithm!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_9ISO_B_V2_rep3.json", b3)

# 10 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7],
                    # :Be9, :B11, :F19, :Ne20, :Na23,
                    # :Mg24, :Al27, :Si28, :P31, :S32,
                    # :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    # :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    # :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    # :Ge74, :As75, :Se80, :Br79, :Kr84,
                    # :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    # :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

a1 = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
end
BenchmarkTools.save("BM1_Kipp_10ISO_A_V2.json", a1)

# a2 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_10ISO_A_V2_rep2.json", a2)

# a3 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_10ISO_A_V2_rep3.json", a3)


b1 = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end
BenchmarkTools.save("BM1_Kipp_10ISO_B_V2.json", b1)

# b2 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
#     Evolution.thomas_algorithm!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_10ISO_B_V2_rep2.json", b2)

# b3 = @benchmark begin
#     StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
#     Evolution.eval_jacobian_eqs!($sm)
#     Evolution.thomas_algorithm!($sm)
# end
# BenchmarkTools.save("BM1_Kipp_10ISO_B_V2_rep3.json", b3)

##

# 11 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9], # :B11, :F19, :Ne20, :Na23],
                    # :Mg24, :Al27, :Si28, :P31, :S32,
                    # :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    # :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    # :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    # :Ge74, :As75, :Se80, :Br79, :Kr84,
                    # :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    # :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_11ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_11ISO_B_V2.json", b)



# 12 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9, :B11], # :F19, :Ne20, :Na23],
                    # :Mg24, :Al27, :Si28, :P31, :S32,
                    # :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    # :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    # :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    # :Ge74, :As75, :Se80, :Br79, :Kr84,
                    # :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    # :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_12ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_12ISO_B_V2.json", b)


# 13 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9, :B11, :F19], # :Ne20, :Na23],
                    # :Mg24, :Al27, :Si28, :P31, :S32,
                    # :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    # :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    # :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    # :Ge74, :As75, :Se80, :Br79, :Kr84,
                    # :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    # :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_13ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_13ISO_B_V2.json", b)


# 14 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9, :B11, :F19, :Ne20], # :Na23,
                    # :Mg24, :Al27, :Si28, :P31, :S32,
                    # :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    # :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    # :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    # :Ge74, :As75, :Se80, :Br79, :Kr84,
                    # :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    # :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_14ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_14ISO_B_V2.json", b)

##


# 15 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9, :B11, :F19, :Ne20, :Na23],
                    # :Mg24, :Al27, :Si28, :P31, :S32,
                    # :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    # :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    # :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    # :Ge74, :As75, :Se80, :Br79, :Kr84,
                    # :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    # :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_15ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_15ISO_B_V2.json", b)



# 20 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9, :B11, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32],
                    # :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    # :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    # :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    # :Ge74, :As75, :Se80, :Br79, :Kr84,
                    # :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    # :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_20ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_20ISO_B_V2.json", b)



# 25 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9, :B11, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45],
                    # :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    # :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    # :Ge74, :As75, :Se80, :Br79, :Kr84,
                    # :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    # :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_25ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_25ISO_B_V2.json", b)



# 30 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9, :B11, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    :Ti48, :V51, :Cr52, :Mn55, :Fe56],
                    # :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    # :Ge74, :As75, :Se80, :Br79, :Kr84,
                    # :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    # :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_30ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_30ISO_B_V2.json", b)



# 35 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9, :B11, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    :Co59, :Ni58, :Cu63, :Zn64, :Ga69],
                    # :Ge74, :As75, :Se80, :Br79, :Kr84,
                    # :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    # :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_35ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_35ISO_B_V2.json", b)




# 40 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9, :B11, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    :Ge74, :As75, :Se80, :Br79, :Kr84],
                    # :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    # :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_40ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_40ISO_B_V2.json", b)



# 45 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9, :B11, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    :Ge74, :As75, :Se80, :Br79, :Kr84,
                    :Rb85, :Sr88, :Y89, :Zr90, :Nb93],
                    # :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_45ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_45ISO_B_V2.json", b)



# 50 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9, :B11, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    :Ge74, :As75, :Se80, :Br79, :Kr84,
                    :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    :Mo98, :Ru102, :Rh103, :Pd106, :Ag107],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_50ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_50ISO_B_V2.json", b)



# 55 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9, :B11, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    :Ge74, :As75, :Se80, :Br79, :Kr84,
                    :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    :Mo98, :Ru102, :Rh103, :Pd106, :Ag107,
                    :Cd114, :In115, :Sn120, :Sb121, :Te130],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_55ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_55ISO_B_V2.json", b)



# 60 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9, :B11, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    :Ge74, :As75, :Se80, :Br79, :Kr84,
                    :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    :Mo98, :Ru102, :Rh103, :Pd106, :Ag107,
                    :Cd114, :In115, :Sn120, :Sb121, :Te130,
                    :I127, :Xe132, :Cs133, :Ba138, :La139],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_60ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_60ISO_B_V2.json", b)



# 65 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9, :B11, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    :Ge74, :As75, :Se80, :Br79, :Kr84,
                    :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    :Mo98, :Ru102, :Rh103, :Pd106, :Ag107,
                    :Cd114, :In115, :Sn120, :Sb121, :Te130,
                    :I127, :Xe132, :Cs133, :Ba138, :La139,
                    :Ce140, :Pr141, :Nd142, :Sm152, :Eu153],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_65ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_65ISO_B_V2.json", b)



# 70 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9, :B11, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    :Ge74, :As75, :Se80, :Br79, :Kr84,
                    :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    :Mo98, :Ru102, :Rh103, :Pd106, :Ag107,
                    :Cd114, :In115, :Sn120, :Sb121, :Te130,
                    :I127, :Xe132, :Cs133, :Ba138, :La139,
                    :Ce140, :Pr141, :Nd142, :Sm152, :Eu153,
                    :Gd158, :Tb159, :Dy164, :Ho165, :Er166],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_70ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_70ISO_B_V2.json", b)




# 75 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9, :B11, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    :Ge74, :As75, :Se80, :Br79, :Kr84,
                    :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    :Mo98, :Ru102, :Rh103, :Pd106, :Ag107,
                    :Cd114, :In115, :Sn120, :Sb121, :Te130,
                    :I127, :Xe132, :Cs133, :Ba138, :La139,
                    :Ce140, :Pr141, :Nd142, :Sm152, :Eu153,
                    :Gd158, :Tb159, :Dy164, :Ho165, :Er166,
                    :Tm169, :Yb174, :Lu175, :Hf180, :Ta181],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_75ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_75ISO_B_V2.json", b)



# 80 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9, :B11, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    :Ge74, :As75, :Se80, :Br79, :Kr84,
                    :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    :Mo98, :Ru102, :Rh103, :Pd106, :Ag107,
                    :Cd114, :In115, :Sn120, :Sb121, :Te130,
                    :I127, :Xe132, :Cs133, :Ba138, :La139,
                    :Ce140, :Pr141, :Nd142, :Sm152, :Eu153,
                    :Gd158, :Tb159, :Dy164, :Ho165, :Er166,
                    :Tm169, :Yb174, :Lu175, :Hf180, :Ta181,
                    :W184, :Re187, :Os192, :Ir193, :Pt195],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_80ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_80ISO_B_V2.json", b)



# 85 ISOTOPES

varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1, :He4, :C12, :N14, :O16,
                    :D2, :T3, :He3, :Li3, :Li7,
                    :Be9, :B11, :F19, :Ne20, :Na23,
                    :Mg24, :Al27, :Si28, :P31, :S32,
                    :Cl35, :Ar36, :K39, :Ca40, :Sc45,
                    :Ti48, :V51, :Cr52, :Mn55, :Fe56,
                    :Co59, :Ni58, :Cu63, :Zn64, :Ga69,
                    :Ge74, :As75, :Se80, :Br79, :Kr84,
                    :Rb85, :Sr88, :Y89, :Zr90, :Nb93,
                    :Mo98, :Ru102, :Rh103, :Pd106, :Ag107,
                    :Cd114, :In115, :Sn120, :Sb121, :Te130,
                    :I127, :Xe132, :Cs133, :Ba138, :La139,
                    :Ce140, :Pr141, :Nd142, :Sm152, :Eu153,
                    :Gd158, :Tb159, :Dy164, :Ho165, :Er166,
                    :Tm169, :Yb174, :Lu175, :Hf180, :Ta181,
                    :W184, :Re187, :Os192, :Ir193, :Pt195,
                    :Au197, :Hg202, :Tl205, :Pb208, :Bi209],
                    [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
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

BenchmarkTools.save("BM1_Kipp_85ISO_A_V2.json", a)

b = @benchmark begin
    StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    Evolution.eval_jacobian_eqs!($sm)
    Evolution.thomas_algorithm!($sm)
end

BenchmarkTools.save("BM1_Kipp_85ISO_B_V2.json", b)









