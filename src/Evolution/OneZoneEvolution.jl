"""
    get_dt_next(oz)

Compute the next one-zone timestep from the change in hydrogen abundance and
the configured timestep growth limits.
"""
function get_dt_next(oz::OneZone)
    dt_next = oz.props.dt
    X = get_value(oz.props.xa[oz.network.xa_index[:H1]])
    Xold = get_value(oz.prv_step_props.xa[oz.network.xa_index[:H1]])
    dt_nextX = dt_next * oz.opt.timestep.delta_Xc_limit / abs(X - Xold)
    minimum_dt = dt_next * oz.opt.timestep.dt_max_decrease
    proposed_dt = min(oz.opt.timestep.dt_max_increase * dt_next, dt_nextX)
    return max(proposed_dt, minimum_dt)
end

"""
    prepare_step!(oz::OneZone, evolver)

Restore the last accepted one-zone model and initialize the timestep for the
next nonlinear trial.
"""
function prepare_step!(oz::OneZone, evolver::DefaultEvolver)
    if evolver.starting_new_step
        evolver.base_timestep = oz.props.dt_next
        evolver.timestep_retry_count = 0
        evolver.starting_new_step = false
    end

    StellarModels.cycle_props!(oz)
    oz.props = deepcopy(oz.prv_step_props)
    oz.props.dt = oz.prv_step_props.dt_next
end

"""
    accept_step!(oz::OneZone, evolver, plotter)

Accept a converged one-zone trial, advance its age and model number, and update
its output and plots.
"""
function accept_step!(oz::OneZone, evolver::DefaultEvolver, plotter)
    evolver.timestep_retry_count = 0
    evolver.starting_new_step = true
    oz.props.time += oz.props.dt
    oz.props.model_number += 1

    StellarModels.evaluate_one_zone_properties!(oz, oz.props)
    StellarModels.write_data(oz)
    StellarModels.write_terminal_info(oz)
    update_plotter!(plotter, oz)
end

"""
    retry_step!(oz::OneZone, evolver, result)

Restore a failed one-zone trial and reduce its timestep. Returns `true` while
another retry is available and `false` once retries are exhausted.
"""
function retry_step!(oz::OneZone, evolver::DefaultEvolver, ::NamedTuple)
    StellarModels.uncycle_props!(oz)

    if evolver.timestep_retry_count < evolver.max_timestep_retries
        evolver.timestep_retry_count += 1
        oz.props.dt_next *= oz.opt.timestep.dt_retry_decrease
        println("Retrying timestep " *
                "$(evolver.timestep_retry_count)/$(evolver.max_timestep_retries)")
        return true
    end

    return false
end

"""
    do_one_zone_burn!(oz; plotter=Plotting.NullPlotter(),
                      max_retries_in_a_row=10,
                      evolver=DefaultEvolver())

Evolve the one-zone burn until a termination condition is reached or the
nonlinear solver exhausts its timestep retries.
"""
function do_one_zone_burn!(oz::OneZone;
                           plotter::TPLOTTER=Plotting.NullPlotter(),
                           max_retries_in_a_row=10,
                           evolver=DefaultEvolver(
                               max_timestep_retries=max_retries_in_a_row,
                           )) where {TPLOTTER<:AbstractPlotter}
    evolver.timestep_retry_count = 0
    evolver.base_timestep = 0.0
    evolver.use_remeshing = false
    evolver.starting_new_step = true
    evolver.skip_remesh_steps_remaining = 0
    StellarModels.create_output_files!(oz)
    StellarModels.evaluate_one_zone_properties!(oz, oz.props)

    trial_active = false
    try
        while true
            trial_active = true
            prepare_step!(oz, evolver)
            result = solve_step!(evolver.solver, oz, evolver)

            if !result.converged
                retrying = retry_step!(oz, evolver, result)
                trial_active = false
                retrying && continue
                println("Terminating evolution: ", result.reason)
                break
            end

            trial_active = false
            accept_step!(oz, evolver, plotter)

            if oz.props.model_number > oz.opt.termination.max_model_number
                StellarModels.write_terminal_info(oz; now=true)
                println("Reached maximum model number")
                break
            end
            if oz.props.T > oz.opt.termination.max_center_T
                StellarModels.write_terminal_info(oz; now=true)
                println("Reached maximum temperature")
                break
            end
            if get_value(oz.props.xa[oz.network.xa_index[:H1]]) <
               oz.opt.termination.min_center_X
                StellarModels.write_terminal_info(oz; now=true)
                println("Reached minimum hydrogen abundance")
                break
            end

            hook_reason = hook_termination_condition(oz, evolver)
            if hook_reason !== nothing
                StellarModels.write_terminal_info(oz; now=true)
                println("Terminating evolution: ", hook_reason)
                break
            end

            oz.props.dt_next = get_dt_next(oz)
            hook_dt_next = hook_get_dt_next(oz, evolver)
            hook_dt_next !== nothing && (oz.props.dt_next = hook_dt_next)
        end
    finally
        trial_active && StellarModels.uncycle_props!(oz)
        StellarModels.shut_down_IO!(oz)
    end
    return nothing
end
