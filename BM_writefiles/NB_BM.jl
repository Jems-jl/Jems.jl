# Import modules

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


x_start = 5
x_end   = 60
n_data  = 20

x_flt = 10 .^(range(log10(x_start), stop = log10(x_end), length = n_data));      # x_flt = float values of normal x-values
x_rnd = round.(Int, x_flt);                                                      # x_rnd = rounded x-values in normal scale
x_log = log10.(x_rnd);                                                           # x_log = x-values in log scale


Isotopes = [:H1, :He4, :C12, :N14, :O16,
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
:Au197, :Hg202, :Tl205, :Pb208, :Bi209,
:Th232, :U238, :Np237, :Pu244, :O15,
:B15, :B16, :C11, :C13, :C14,
:Li4, :Li5, :F18, :Ne32, :Mg26,
:Lv290, :Lv291, :Lv292, :Lv293, :Pb206,
:Pb207, :Bi207, :Bi208, :Po208, :Po209,
:Po210, :Rn210, :Rn211, :Rn222, :Ra226,
:Ra228, :U233, :U234, :U235, :U236]

# Model creation

for n in x_rnd

    println("Amount of isotopes: ", n)

    varnames = [:lnρ, :lnT, :lnr, :lum]
    varscaling = [:log, :log, :log, :maxval]
    structure_equations = [Evolution.equationHSE, Evolution.equationT,
                           Evolution.equationContinuity, Evolution.equationLuminosity]
    remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                              StellarModels.split_lnT, StellarModels.split_xa]
    net = NuclearNetwork(Isotopes[1:n], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
    nz = 1000
    nextra = 100
    eos = EOS.IdealEOS(true)
    opacity = Opacity.SimpleElectronScatteringOpacity()
    turbulence = Turbulence.BasicMLT(1.0)
    sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                        nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);

    n = 3
    StellarModels.n_polytrope_initial_condition!(n, sm, nz, 0.7154, 0.0142, 0.0003, Chem.abundance_lists[:ASG_09], MSUN,
                                                 100 * RSUN; initial_dt=10 * SECYEAR)
    StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
    StellarModels.cycle_props!(sm);
    StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
    StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)  # or do StellarModels.remesher!(sm);
    StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
    StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
    StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

    println("Start benchmarks")

    # ----------

    folder_path_eval_and_jac = "/Users/evakuipers/Jems.jl/BM_outputfolders/BM_NB_eval_and_jac_output" 

    a = @benchmark begin
        StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
        Evolution.eval_jacobian_eqs!($sm)
    end

    file_path_eval_and_jac = joinpath(folder_path_eval_and_jac, "BM_NB_eval_and_jac_ISO$(n).json")
    BenchmarkTools.save(file_path_eval_and_jac, a)

    println("Benchmark 1 done")

    # ----------

    folder_path_eval = "/Users/evakuipers/Jems.jl/BM_outputfolders/BM_NB_eval_output"  

    b = @benchmark begin
        StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
    end

    file_path_eval = joinpath(folder_path_eval, "BM_NB_eval_ISO$(n).json")
    BenchmarkTools.save(file_path_eval, b)

    println("Benchmark 2 done")

    # ----------

    folder_path_jac = "/Users/evakuipers/Jems.jl/BM_outputfolders/BM_NB_jac_output"  

    c = @benchmark begin
        setup = (StellarModels.evaluate_stellar_model_properties!($sm, $sm.props))
        Evolution.eval_jacobian_eqs!($sm)
    end

    file_path_jac = joinpath(folder_path_jac, "BM_NB_jac_ISO$(n).json")
    BenchmarkTools.save(file_path_jac, c)

    println("Benchmark 3 done")

    # ----------

    folder_path_solver = "/Users/evakuipers/Jems.jl/BM_outputfolders/BM_NB_solver_output"  

    d = @benchmark begin
        setup = (StellarModels.evaluate_stellar_model_properties!($sm, $sm.props), Evolution.eval_jacobian_eqs!($sm))
        Evolution.thomas_algorithm!($sm)
    end

    file_path_solver = joinpath(folder_path_solver, "BM_NB_solver_ISO$(n).json")
    BenchmarkTools.save(file_path_solver, d)

    println("Benchmark 4 done")

    # ----------

    folder_path_eval_jac_solver = "/Users/evakuipers/Jems.jl/BM_outputfolders/BM_NB_eval_jac_solver_output"  

    e = @benchmark begin
        StellarModels.evaluate_stellar_model_properties!($sm, $sm.props)
        Evolution.eval_jacobian_eqs!($sm)
        Evolution.thomas_algorithm!($sm)
    end

    file_path_eval_jac_solver = joinpath(folder_path_eval_jac_solver, "BM_NB_eval_jac_solver_ISO$(n).json")
    BenchmarkTools.save(file_path_eval_jac_solver, e)

    println("Benchmark 5 done")

end

