Base.@kwdef struct TimestepRetryPolicy <: AbstractRetryPolicy
    maximum_retries::Int = 10
    reduction_factor::Float64 = 0.5
    retry_without_remeshing::Bool = true
end

function retry_step!(retry_policy::TimestepRetryPolicy, sm::StellarModel, evolver::StellarEvolver, state::StellarEvolutionState, result::StepResult)
    StellarModels.uncycle_props!(sm)  # reset props to what prv_step_props contains, ie mimic state at end of previous step
    
    # Try another timestep using the current mesh mode
    if state.timestep_retry_count < retry_policy.max_timestep_retries
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
