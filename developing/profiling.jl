using BenchmarkTools
using Jems.Chem
using Jems.Constants
using Jems.EOS
using Jems.Opacity
using Jems.NuclearNetworks
using Jems.StellarModels
using Jems.Evolution
using Jems.ReactionRates
using Jems.Turbulence
using Profile
using PProf
using Jems.DualSupport
using ForwardDiff

##
varnames = [:lnρ, :lnT, :lnr, :lum]
varscaling = [:log, :log, :log, :maxval]
structure_equations = [Evolution.equationHSE, Evolution.equationT,
                       Evolution.equationContinuity, Evolution.equationLuminosity]
remesh_split_functions = [StellarModels.split_lnr_lnρ, StellarModels.split_lum,
                          StellarModels.split_lnT, StellarModels.split_xa]
net = NuclearNetwork([:H1,:He4,:C12, :N14, :O16], [(:kipp_rates, :kipp_pp), (:kipp_rates, :kipp_cno)])
nz = 1000
nextra = 100
eos = EOS.IdealEOS(false)
opacity = Opacity.SimpleElectronScatteringOpacity()
turbulence = Turbulence.BasicMLT(1.0)
sm = StellarModel(varnames, varscaling, structure_equations, Evolution.equation_composition,
                  nz, nextra, remesh_split_functions, net, eos, opacity, turbulence);

##
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
##
function copy_partials!(m::AbstractModel, k::Int)
    jacobian_Lk = m.solver_data.jacobian_L[k]
    jacobian_Dk = m.solver_data.jacobian_D[k]
    jacobian_Uk = m.solver_data.jacobian_U[k]
    eqs_duals = m.solver_data.eqs_duals
    for i = 1:(m.nvars)
        for j = 1:(m.nvars)
            jacobian_Lk[i, j] = eqs_duals[k, i].partials[j]
            jacobian_Dk[i, j] = eqs_duals[k, i].partials[j + m.nvars]
            jacobian_Uk[i, j] = eqs_duals[k, i].partials[j + 2 * m.nvars]
        end
        # populate the eqs_numbers with relevant entries (will be RHS for linear solver)
        m.solver_data.eqs_numbers[(k - 1) * m.nvars + i] = eqs_duals[k, i].value
    end
end
##
function copy_partials!(m::AbstractModel, k::Int)
    jacobian_Lk = m.solver_data.jacobian_L[k]
    jacobian_Dk = m.solver_data.jacobian_D[k]
    jacobian_Uk = m.solver_data.jacobian_U[k]
    eqs_duals = m.solver_data.eqs_duals
    for i = 1:(m.nvars)
        block_select = eqs_duals[k, i]
        for j = 1:(m.nvars)
            jacobian_Lk[i, j] = block_select.partials[j]
            jacobian_Dk[i, j] = block_select.partials[j + m.nvars]
            jacobian_Uk[i, j] = block_select.partials[j + 2 * m.nvars]
        end
        # populate the eqs_numbers with relevant entries (will be RHS for linear solver)
        m.solver_data.eqs_numbers[(k - 1) * m.nvars + i] = eqs_duals[k, i].value
    end
end

##
Profile.clear()
Profile.@profile begin
    # StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
    for i = 1:10000
        copy_partials!(sm, 1)
    end
    # Evolution.thomas_algorithm!(sm)
end
PProf.pprof()

##
@benchmark Evolution.eval_jacobian_eqs_row!($sm, 2)

##
@benchmark Evolution.eval_cell_eqs!($sm, 2)

##
@benchmark copy_partials!($sm, 2)

##
Profile.clear()
Profile.@profile begin
    for i = 1:10000
        my_eval_jacobian_eqs_row!(sm, 2)
    end
end
PProf.pprof()

##
@benchmark Evolution.eval_jacobian_eqs!($sm)
