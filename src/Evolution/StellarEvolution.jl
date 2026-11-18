"""
    get_dt_next(sm)

Compute the next stellar timestep from changes in surface radius, temperature,
and central hydrogen abundance, subject to the configured timestep growth and
reduction limits.
"""
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

"""
    compute_starting_model_properties!(sm, evolver=DefaultEvolver())

Evaluate the initial stellar model and populate the property buffers used to
start the first evolutionary step.
"""
function compute_starting_model_properties!(sm::StellarModel,
                                            evolver::DefaultEvolver=DefaultEvolver())
    hook_before_evaluate_stellar_model_properties!(sm, sm.props, evolver)
    StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
    hook_after_evaluate_stellar_model_properties!(sm, sm.props, evolver)
    StellarModels.cycle_props!(sm)
    StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
    StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)
    hook_before_evaluate_stellar_model_properties!(sm, sm.start_step_props, evolver)
    StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
    hook_after_evaluate_stellar_model_properties!(sm, sm.start_step_props, evolver)
    StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
    StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)
end

"""
    prepare_step!(model, evolver)

Prepare a trial step from the last accepted model. For a `StellarModel`, this
also selects the remeshed or original-mesh retry path.
"""
function prepare_step!(sm::StellarModel, evolver::DefaultEvolver)
    if evolver.starting_new_step
        evolver.base_timestep = sm.props.dt_next
        if evolver.skip_remesh_steps_remaining > 0
            evolver.use_remeshing = false
            evolver.skip_remesh_steps_remaining -= 1
        else
            evolver.use_remeshing = sm.opt.remesh.do_remesh
        end
        evolver.timestep_retry_count = 0
        evolver.starting_new_step = false
    end
    StellarModels.cycle_props!(sm)

    if sm.opt.remesh.do_remesh && evolver.use_remeshing
        StellarModels.remesher!(sm)
    else
        StellarModels.copy_mesh_properties!(sm, sm.start_step_props,
                                            sm.prv_step_props)
    end

    StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
    hook_before_evaluate_stellar_model_properties!(sm, sm.start_step_props, evolver)
    StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
    hook_after_evaluate_stellar_model_properties!(sm, sm.start_step_props, evolver)
    sm.start_step_props.dt = sm.prv_step_props.dt_next
    StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
    StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)
end

"""
    accept_step!(model, evolver, plotter)

Accept a converged trial, advance its age and model number, write output, and
update `plotter`.
"""
function accept_step!(sm::StellarModel, evolver::DefaultEvolver, plotter)
    evolver.timestep_retry_count = 0
    evolver.starting_new_step = true
    sm.props.time += sm.props.dt
    sm.props.model_number += 1

    hook_before_evaluate_stellar_model_properties!(sm, sm.props, evolver)
    StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
    hook_after_evaluate_stellar_model_properties!(sm, sm.props, evolver)
    StellarModels.write_data(sm, eltype(sm.props.ind_vars))
    StellarModels.write_terminal_info(sm)
    update_plotter!(plotter, sm)
end

"""
    retry_step!(model, evolver, result)

Restore a failed trial and select the next recovery attempt. Returns `true`
when another attempt should be made and `false` when retries are exhausted.
"""
function retry_step!(sm::StellarModel, evolver::DefaultEvolver, ::NamedTuple)
    StellarModels.uncycle_props!(sm)

    if evolver.timestep_retry_count < evolver.max_timestep_retries
        evolver.timestep_retry_count += 1
        sm.props.dt_next *= sm.opt.timestep.dt_retry_decrease
        mesh = evolver.use_remeshing ? "remeshed" : "original"
        println("Retrying $mesh mesh with timestep retry " *
                "$(evolver.timestep_retry_count)/$(evolver.max_timestep_retries)")
        return true
    end

    if evolver.use_remeshing
        evolver.use_remeshing = false
        evolver.timestep_retry_count = 0
        evolver.skip_remesh_steps_remaining = evolver.skip_remesh_steps
        sm.props.dt_next = evolver.base_timestep
        println("Remeshed attempts exhausted; retrying the original mesh " *
                "with the base timestep")
        if evolver.skip_remesh_steps > 0
            println("Remeshing will remain disabled for the next " *
                    "$(evolver.skip_remesh_steps) accepted steps")
        end
        return true
    end

    return false
end



"""
    do_evolution_loop!(sm; plotter=Plotting.NullPlotter(),
                       max_retries_in_a_row=10,
                       evolver=DefaultEvolver())

Evolve `sm` until a built-in or user-defined termination condition is reached,
or until the nonlinear solver exhausts its recovery attempts.
"""
function do_evolution_loop!(sm::StellarModel;
                            plotter::TPLOTTER=Plotting.NullPlotter(),
                            max_retries_in_a_row=10,
                            evolver=DefaultEvolver(
                                max_timestep_retries=max_retries_in_a_row,
                            )) where {TPLOTTER<:AbstractPlotter}
    evolver.timestep_retry_count = 0
    evolver.base_timestep = 0.0
    evolver.use_remeshing = sm.opt.remesh.do_remesh
    evolver.starting_new_step = true
    evolver.skip_remesh_steps_remaining = 0
    StellarModels.create_output_files!(sm, eltype(sm.props.ind_vars))
    compute_starting_model_properties!(sm, evolver)

    trial_active = false
    try
        while true
            trial_active = true
            prepare_step!(sm, evolver)
            result = solve_step!(evolver.solver, sm, evolver)

            if !result.converged
                retrying = retry_step!(sm, evolver, result)
                trial_active = false
                retrying && continue
                println("Terminating evolution: ", result.reason)
                break
            end

            trial_active = false
            accept_step!(sm, evolver, plotter)

            # Check termination conditions.
            if sm.props.model_number > sm.opt.termination.max_model_number
                StellarModels.write_terminal_info(sm; now=true)
                println("Reached maximum model number")
                break
            end
            if exp(get_value(sm.props.lnT[1])) >
               sm.opt.termination.max_center_T
                StellarModels.write_terminal_info(sm; now=true)
                println("Reached maximum central temperature")
                break
            end
            if get_value(sm.props.xa[1, sm.network.xa_index[:H1]]) <
               sm.opt.termination.min_center_X
                StellarModels.write_terminal_info(sm; now=true)
                println("Reached minimum central hydrogen abundance")
                break
            end

            # User-provided termination hooks run after the regular checks.
            hook_reason = hook_termination_condition(sm, evolver)

            if hook_reason !== nothing
                StellarModels.write_terminal_info(sm; now=true)
                println("Terminating evolution: ", hook_reason)
                break
            end

            sm.props.dt_next = get_dt_next(sm)
            hook_dt_next = hook_get_dt_next(sm, evolver)
            hook_dt_next !== nothing && (sm.props.dt_next = hook_dt_next)
        end
    finally
        trial_active && StellarModels.uncycle_props!(sm)
        StellarModels.shut_down_IO!(sm)
    end
    return nothing
end
