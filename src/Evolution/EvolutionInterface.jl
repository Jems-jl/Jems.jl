"""
    AbstractEvolver

Abstract struct for objects that configure an evolution of the different StellarModels.
"""
abstract type AbstractEvolver end

"""
    AbstractNonlinearSolver

Abstract struct for nonlinear solver strategies used to solve one evolutionary step.
"""
abstract type AbstractNonlinearSolver end

"""
    DefaultEvolver

Default evolution controller shared currently by stellar and one-zone models. It stores
the nonlinear solver, the permitted number of timestep retries, and the
mutable state needed while recovering a failed step.

Construct one with `DefaultEvolver(; solver=NewtonSolver(),
max_timestep_retries=10, skip_remesh_steps=0)`.
"""
mutable struct DefaultEvolver{TSOLVER<:AbstractNonlinearSolver} <: AbstractEvolver
    solver::TSOLVER
    max_timestep_retries::Int
    timestep_retry_count::Int
    base_timestep::Float64
    use_remeshing::Bool
    starting_new_step::Bool
    skip_remesh_steps::Int
    skip_remesh_steps_remaining::Int
end

"""
    DefaultEvolver(; solver=NewtonSolver(), max_timestep_retries=10,
                   skip_remesh_steps=0)

Construct the default evolver using `solver` for each nonlinear step. After a
failed remeshed step falls back to the original mesh, `skip_remesh_steps`
controls how many subsequent accepted steps also keep the original mesh.
"""
function DefaultEvolver(; solver=NewtonSolver(), max_timestep_retries=10,
                        skip_remesh_steps=0)
    skip_remesh_steps >= 0 || throw(ArgumentError("skip_remesh_steps must be nonnegative"))
    return DefaultEvolver(solver, max_timestep_retries, 0, 0.0, true, true,
                          skip_remesh_steps, 0)
end

# Hook functions 
# These functions in general return nothing but user can modify them from outside to
# add extra functionalities including extra termination conditions, extra time-step controls, etc 

"""
    hook_termination_condition(model, evolver)

Optional user-defined termination check. Return `nothing` to continue evolving,
or return a value describing why the evolution should stop.
"""
hook_termination_condition(model, evolver::AbstractEvolver) = nothing

"""
    hook_before_evaluate_stellar_model_properties!(model, props, evolver)

Optional hook run immediately before `props` are evaluated. The default method
does nothing.
"""
hook_before_evaluate_stellar_model_properties!(model, props,
                                          evolver::AbstractEvolver) = nothing

"""
    hook_after_evaluate_stellar_model_properties!(model, props, evolver)

Optional hook run after `props` have been evaluated successfully. The default
method does nothing.
"""
hook_after_evaluate_stellar_model_properties!(model, props,
                                         evolver::AbstractEvolver) = nothing

"""
    hook_get_dt_next(model, evolver)

Return `nothing` to use the default timestep controller, or return the desired
next timestep using user defined controls.
"""
hook_get_dt_next(model, evolver::AbstractEvolver) = nothing
