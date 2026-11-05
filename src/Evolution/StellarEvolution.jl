"""Compute the timestep proposed for the next stellar model."""
function get_dt_next(sm::StellarModel)
    dt_next = sm.props.dt

    radius = exp(get_value(sm.props.lnr[sm.props.nz]))
    old_radius = exp(get_value(sm.prv_step_props.lnr[sm.prv_step_props.nz]))
    delta_radius = abs(radius - old_radius) / radius

    temperature = exp(get_value(sm.props.lnT[sm.props.nz]))
    old_temperature = exp(get_value(sm.prv_step_props.lnT[sm.prv_step_props.nz]))
    delta_temperature = abs(temperature - old_temperature) / temperature

    hydrogen = get_value(sm.props.xa[1, sm.network.xa_index[:H1]])
    old_hydrogen = get_value(sm.prv_step_props.xa[1, sm.network.xa_index[:H1]])
    delta_hydrogen = abs(hydrogen - old_hydrogen)

    dt_radius = dt_next * sm.opt.timestep.delta_R_limit / delta_radius
    dt_temperature = dt_next * sm.opt.timestep.delta_Tc_limit / delta_temperature
    dt_hydrogen = dt_next * sm.opt.timestep.delta_Xc_limit / delta_hydrogen
    minimum_dt = dt_next * sm.opt.timestep.dt_max_decrease
    proposed_dt = min(sm.opt.timestep.dt_max_increase * dt_next,
                      dt_radius, dt_temperature, dt_hydrogen)
    return max(proposed_dt, minimum_dt)
end

function compute_starting_model_properties!(sm::StellarModel)
    StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
    StellarModels.cycle_props!(sm)
    StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
    StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)
    StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
    StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
    StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)
end

function initialize_base_step!(sm::StellarModel, state::StellarEvolutionState)
    state.base_timestep = sm.props.dt_next
    state.use_remeshing = sm.opt.remesh.do_remesh
    state.timestep_retry_count = 0
    state.starting_new_step = false
end

function reset_recovery_after_success!(state::StellarEvolutionState)
    state.timestep_retry_count = 0
    state.starting_new_step = true
end

function prepare_step!(sm::StellarModel, ::StellarEvolver,
                       state::StellarEvolutionState)
    state.starting_new_step && initialize_base_step!(sm, state)
    StellarModels.cycle_props!(sm)

    if sm.opt.remesh.do_remesh && state.use_remeshing
        StellarModels.remesher!(sm)
    else
        StellarModels.copy_mesh_properties!(sm, sm.start_step_props,
                                            sm.prv_step_props)
    end

    StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
    StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
    sm.start_step_props.dt = sm.prv_step_props.dt_next
    StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
    StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)
end

function accept_step!(sm::StellarModel, ::StellarEvolver,
                      state::StellarEvolutionState, plotter)
    reset_recovery_after_success!(state)
    sm.props.time += sm.props.dt
    sm.props.model_number += 1

    StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
    StellarModels.write_data(sm, eltype(sm.props.ind_vars))
    StellarModels.write_terminal_info(sm)
    update_plotter!(plotter, sm)
end

"""Restore a failed trial and select the next timestep/mesh attempt."""
function retry_step!(sm::StellarModel, evolver::StellarEvolver,
                     state::StellarEvolutionState, ::StepResult)
    StellarModels.uncycle_props!(sm)

    if state.timestep_retry_count < evolver.max_timestep_retries
        state.timestep_retry_count += 1
        sm.props.dt_next *= sm.opt.timestep.dt_retry_decrease
        mesh = state.use_remeshing ? "remeshed" : "original"
        println("Retrying $mesh mesh with timestep retry " *
                "$(state.timestep_retry_count)/$(evolver.max_timestep_retries)")
        return true
    end

    if state.use_remeshing
        state.use_remeshing = false
        state.timestep_retry_count = 0
        sm.props.dt_next = state.base_timestep
        println("Remeshed attempts exhausted; retrying the original mesh " *
                "with the base timestep")
        return true
    end

    return false
end

function termination_reason(sm::StellarModel, evolver::StellarEvolver,
                            state::StellarEvolutionState)
    sm.props.model_number > sm.opt.termination.max_model_number &&
        return :maximum_model_number
    exp(get_value(sm.props.lnT[1])) > sm.opt.termination.max_center_T &&
        return :maximum_center_temperature
    get_value(sm.props.xa[1, sm.network.xa_index[:H1]]) <
        sm.opt.termination.min_center_X && return :minimum_center_hydrogen

    for condition in evolver.extra_termination_conditions
        reason = condition(sm, evolver, state)
        reason === nothing || return reason
    end
    return nothing
end

function evolve!(sm::StellarModel, evolver::StellarEvolver;
                 plotter=Plotting.NullPlotter())
    state = StellarEvolutionState()
    StellarModels.create_output_files!(sm, eltype(sm.props.ind_vars))
    compute_starting_model_properties!(sm)

    trial_active = false
    try
        while true
            trial_active = true
            prepare_step!(sm, evolver, state)
            result = solve_step!(evolver.solver, sm, evolver, state)

            if !result.converged
                retrying = retry_step!(sm, evolver, state, result)
                trial_active = false
                retrying && continue
                println("Terminating evolution: ", result.reason)
                break
            end

            trial_active = false
            accept_step!(sm, evolver, state, plotter)
            reason = termination_reason(sm, evolver, state)
            if reason !== nothing
                StellarModels.write_terminal_info(sm; now=true)
                println("Terminating evolution: ", reason)
                break
            end
            sm.props.dt_next = get_dt_next(sm)
        end
    finally
        trial_active && StellarModels.uncycle_props!(sm)
        StellarModels.shut_down_IO!(sm)
    end
    return nothing
end

"""Compatibility entry point for stellar evolution."""
function do_evolution_loop!(sm::StellarModel;
                            plotter::TPLOTTER=Plotting.NullPlotter(),
                            max_retries_in_a_row=10,
                            evolver=StellarEvolver(
                                max_timestep_retries=max_retries_in_a_row,
                            )) where {TPLOTTER<:AbstractPlotter}
    return evolve!(sm, evolver; plotter=plotter)
end
