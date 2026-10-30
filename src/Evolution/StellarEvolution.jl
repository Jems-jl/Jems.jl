endinclude("EvolutionInterface.jl")
include("NonlinearSolvers.jl")
include("RetryPolicies.jl")
"""
    compute_starting_model_properties(sm::StellarModel)

Computes all stellar model properties based on the current state of the star.
This sets the model up to start the evolution loop, or test the evaluation of
the equations and the linear solver.
"""
function compute_starting_model_properties!(sm::StellarModel)
    StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
    StellarModels.cycle_props!(sm);
    StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
    StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)
    StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
    StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
    StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)
end

function initialize_base_step!(sm::StellarModel, state::EvolutionState)
    state.base_timestep = sm.props.dt_next
    state.use_remeshing = sm.opt.remesh.do_remesh 
    state.timestep_retry_count = 0 
    state.starting_new_step = false
end 
function reset_recovery_after_success!( state::EvolutionState,)
    state.timestep_retry_count = 0
    state.use_remeshing = true
    state.starting_new_step = true
end

function prepare_step!(sm::StellarModel, evolver::AbstractEvolver, state::EvolutionState) #everything related to step preapration 
    if state.starting_new_step
        initialize_base_step!(sm, state)
    end
    StellarModels.cycle_props!(sm) # move props of previous step to prv_step_props of current step

    # checking/setting the remeshing status for the step 
    should_remesh = sm.opt.remesh.do_remesh && state.use_remeshing # set to true if toml says do remesh and no remesh is not forced 
    # either remesh, or copy over from prv_step_props

    remesh_model!(evolver.remesh, sm, state)
    # time derivatives in the equations use the remeshed info, so
    # save start_step_props before we attempt any newton solver
    StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
    StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
    sm.start_step_props.dt = sm.prv_step_props.dt_next  # dt of this step becomes dt_next of previous
    StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
    StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)
end 


function accept_step!(sm::StellarModel, evolver::AbstractEvolver, state::EvolutionState, plotter,)
    reset_recovery_after_success!(state)
    # increment age and model number since we accept the step.
    sm.props.time += sm.props.dt 
    sm.props.model_number += 1 

    # write state in sm.props and potential history/profiles.
    StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
    StellarModels.write_data(sm, eltype(sm.props.ind_vars))
    StellarModels.write_terminal_info(sm)

    update_plotter!(plotter, sm)
end 


function evolve!(sm::StellarModel, evolver::StellarEvolver;plotter=Plotting.NullPlotter(),)
    state = StellarEvolutionState()
    StellarModels.create_output_files!(sm, eltype(sm.props.ind_vars))
    compute_starting_model_properties!(sm)

    timestep_trial_active = false # to check if we are doing timestep trial 
    try 
        while true 
            timestep_trial_active = true
            prepare_step!(sm, evolver,state)
            result = solve_step!(evolver.solver, sm, evolver, state)

            if !result.converged
                retrying = retry_step!(evolver.retry, sm, evolver, state, result)
                timestep_trial_active = false
                if retrying 
                    continue 
                end 

                println("Terminating evolution: ", result.reason)
                break
            end 
            # when converged
            timestep_trial_active = false
            accept_step!(sm, evolver, state, plotter,)
            reason = termination_reason(evolver.termination, sm, evolver, state)

            if reason !== nothing 
                StellarModels.write_terminal_info(sm; now=true)
                println("Terminating evolution: ", reason)
                break
            end 
            sm.props.dt_next = get_dt_next(evolver.timestep, sm, evolver, state)
        end 
    finally 
        if timestep_trial_active
         StellarModels.uncycle_props!(sm)
        end
        StellarModels.shut_down_IO!(sm)
    end 

    return nothing 
end 
"""
    do_evolution_loop(sm::StellarModel)

Performs the main evolutionary loop of the input StellarModel `sm`. It continues taking steps until one of the
termination criteria is reached (defined in `sm.opt.termination`).
"""
function do_evolution_loop!(sm::StellarModel; plotter::TPLOTTER = Plotting.NullPlotter(),
                                max_retries_in_a_row = 10,  evolver=StellarEvolver(max_timestep_retries=max_retries_in_a_row,)) where{TPLOTTER<:AbstractPlotter}
        return evolve!(sm, evolver; plotter=plotter,)
end 