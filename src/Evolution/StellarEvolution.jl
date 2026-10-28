abstract type AbstractEvolver end  #Contains every struct related to evolution 

@kwdef struct StellarEvolver <: AbstractEvolver #describes how to evolove 
    max_timestep_retries::Int = 10 
    # residual_tolerannce
    # correction_tolerance 

end 

Base.@kwdef mutable struct EvolutionState
    timestep_retry_count::Int = 0
    base_timestep::Float64 = 0.0
    use_remeshing::Bool = true
    starting_new_step::Bool = true
end

struct StepResult
    converged::Bool 
    reason::Symbol
    iterations::Int
end 
"""
    get_dt_next(sm::StellarModel)

Computes the timestep of the next evolutionary step to be taken by the StellarModel `sm` by considering all timestep
controls (`sm.opt.timestep`).
"""
function get_dt_next(sm::StellarModel)
    dt_next = sm.props.dt  # this it calculated at end of step, so props.dt is the dt we used to do this step
        
    Rsurf = exp(get_value(sm.props.lnr[sm.props.nz]))
    Rsurf_old = exp(get_value(sm.prv_step_props.lnr[sm.prv_step_props.nz]))
    ΔR_div_R = abs(Rsurf - Rsurf_old) / Rsurf

    Tc = exp(get_value(sm.props.lnT[sm.props.nz]))
    Tc_old = exp(get_value(sm.prv_step_props.lnT[sm.prv_step_props.nz]))
    ΔTc_div_Tc = abs(Tc - Tc_old) / Tc

    X = get_value(sm.props.xa[1, sm.network.xa_index[:H1]])
    Xold = get_value(sm.prv_step_props.xa[1, sm.network.xa_index[:H1]])
    ΔX = abs(X - Xold)

    dt_nextR = dt_next * sm.opt.timestep.delta_R_limit / ΔR_div_R
    dt_nextTc = dt_next * sm.opt.timestep.delta_Tc_limit / ΔTc_div_Tc
    dt_nextX = dt_next * sm.opt.timestep.delta_Xc_limit / ΔX

    min_dt = dt_next * sm.opt.timestep.dt_max_decrease
    dt_next = min(sm.opt.timestep.dt_max_increase * dt_next, dt_nextR, dt_nextTc, dt_nextX)
    dt_next = max(dt_next, min_dt)
    return dt_next
end

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
function remesh_model!(sm::StellarModel, ::StellarEvolver) # we change it later 
    StellarModels.remesher!(sm)
end 

function prepare_step!(sm::StellarModel, evolver::AbstractEvolver, state::EvolutionState) #everything related to step preapration 
    if state.starting_new_step
        initialize_base_step!(sm, state)
    end
    StellarModels.cycle_props!(sm) # move props of previous step to prv_step_props of current step

    # checking/setting the remeshing status for the step 
    should_remesh = sm.opt.remesh.do_remesh && state.use_remeshing # set to true if toml says do remesh and no remesh is not forced 
    # either remesh, or copy over from prv_step_props

    if should_remesh
        remesh_model!(sm, evolver)
    else
        StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)
    end 
    # time derivatives in the equations use the remeshed info, so
    # save start_step_props before we attempt any newton solver
    StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
    StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
    sm.start_step_props.dt = sm.prv_step_props.dt_next  # dt of this step becomes dt_next of previous
    StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
    StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)
end 

function solve_step!(sm::StellarModel, evolver::AbstractEvolver, state::EvolutionState)
    # step loop
    sm.solver_data.newton_iters = 0
    max_steps = sm.opt.solver.newton_max_iter
    if (sm.start_step_props.model_number == 0)
        max_steps = sm.opt.solver.newton_max_iter_first_step
    end
    corr = @view sm.solver_data.solver_corr[1:sm.nvars*sm.props.nz]
    equs = @view sm.solver_data.eqs_numbers[1:sm.nvars*sm.props.nz]
    
    # evaluate the equations for the first step
    try 
        StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
        eval_jacobian_eqs!(sm)  # heavy lifting happens here!]
    catch e
        e isa InterruptException && rethrow()
        println("Error while evaluating initial equations")
        showerror(stdout,e)
        println()
        return StepResult(false, :initial_evaluation_error, 0,)
    end 

    for i = 1:max_steps
      
        try

            block_tridiagonal_solver!(sm, sm.solver_data)  # here as well

            (abs_max_corr, i_corr) = findmax(abs, corr)
            signed_max_corr = corr[i_corr]
            corr_nz = div(i_corr - 1, sm.nvars) + 1
            corr_equ = mod1(i_corr, sm.nvars)
            rel_corr = abs_max_corr/eps(sm.props.ind_vars[i_corr])

            # scale correction
            if sm.props.model_number == 0
                correction_limit = sm.opt.solver.initial_model_scale_max_correction
            else
                correction_limit = sm.opt.solver.scale_max_correction
            end
            correction_multiplier = 1.0
            for j in 1:sm.nvars
                sub_corr = @view sm.solver_data.solver_corr[j:sm.nvars:(sm.nvars*(sm.props.nz-1)+j)]
                max_sub_corr = maximum(abs, sub_corr)
                if sm.var_scaling[j] == :log || sm.var_scaling[j] == :unity
                    correction_multiplier = min(correction_multiplier,
                                                    correction_limit/max_sub_corr)
                elseif sm.var_scaling[j] == :maxval
                    sub_ind_vars = @view sm.props.ind_vars[j:sm.nvars:(sm.nvars*(sm.props.nz-1)+j)]
                    max_var_value = maximum(abs, sub_ind_vars)
                    correction_multiplier = min(correction_multiplier,
                                                    max_var_value*correction_limit/max_sub_corr)
                end
            end

            if correction_multiplier < 1.0
                corr .*= correction_multiplier
            end

            # apply correction!
            sm.props.ind_vars[1:sm.nvars*sm.props.nz] .+= corr[1:sm.nvars*sm.props.nz]
            sm.solver_data.newton_iters = i

            # evaluate the equations after correction and get residuals
            StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
            eval_jacobian_eqs!(sm)  # heavy lifting happens here!

            (max_res, i_res) = findmax(abs, equs)
            res_nz = div(i_res - 1, sm.nvars) + 1
            res_equ =mod1(i_res, sm.nvars)

            #reporting
            if sm.opt.solver.report_solver_progress &&
                i % sm.opt.solver.solver_progress_iter == 0
                @show sm.props.model_number, i, rel_corr, signed_max_corr, corr_nz, corr_equ, max_res, res_nz, res_equ
            end
            #check if tolerances are satisfied
            if rel_corr < sm.opt.solver.relative_correction_tolerance &&
                    max_res < sm.opt.solver.maximum_residual_tolerance
                if sm.props.model_number == 0
                    println("Found first model")
                end
                return StepResult(true,:converged,i)
            end
        catch e
            e isa InterruptException && rethrow()
            println("Error during Newton iteration $i")
            showerror(stdout, e)
            println()
            return StepResult(false, :evaluation_error, i)
        end
    end
    return StepResult(false, :iteration_limit, max_steps)
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

function retry_step!(sm::StellarModel, evolver::StellarEvolver, state::EvolutionState, result::StepResult)
    StellarModels.uncycle_props!(sm)  # reset props to what prv_step_props contains, ie mimic state at end of previous step
    
    # Try another timestep using the current mesh mode
    if state.timestep_retry_count < evolver.max_timestep_retries
        state.timestep_retry_count += 1
        sm.props.dt_next *= sm.opt.timestep.dt_retry_decrease
        mesh_mode = state.use_remeshing ? "remeshed" : "original"
        println("Retrying", mesh_mode, " mesh with timestep retry ", state.timestep_retry_count)
        return true 
    
    end 
    if state.use_remeshing # if it is true, then we set things to false, so that things run with with original mesh 
        state.use_remeshing = false
        state.timestep_retry_count = 0 
        sm.props.dt_next = state.base_timestep
        return true 
    end 
    return false
end 

function termination_reason(sm::StellarModel, evolver::AbstractEvolver, state::EvolutionState)
    if sm.props.model_number > sm.opt.termination.max_model_number
        return :maximum_model_number
    end

    if exp(get_value(sm.props.lnT[1])) >
       sm.opt.termination.max_center_T
        return :maximum_center_temperature
    end

    center_hydrogen = get_value(sm.props.xa[1, sm.network.xa_index[:H1]])

    if center_hydrogen < sm.opt.termination.min_center_X
        return :minimum_center_hydrogen
    end

    return nothing
end 

function evolve!(sm::StellarModel, evolver::AbstractEvolver;plotter=Plotting.NullPlotter(),)
    state = EvolutionState()
    StellarModels.create_output_files!(sm, eltype(sm.props.ind_vars))
    compute_starting_model_properties!(sm)

    timestep_trial_active = false # to check if we are doing timestep trial 
    try 
        while true 
            timestep_trial_active = true
            prepare_step!(sm, evolver,state)
            result = solve_step!(sm, evolver, state)

            if !result.converged
                retrying = retry_step!(sm,evolver, state, result)
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
            reason = termination_reason(sm, evolver, state)

            if reason !== nothing 
                StellarModels.write_terminal_info(sm; now=true)
                println("Terminating evolution: ", reason)
                break
            end 
            sm.props.dt_next = get_dt_next(sm)
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