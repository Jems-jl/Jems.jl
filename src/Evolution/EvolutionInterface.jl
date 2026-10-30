abstract type AbstractEvolver end
abstract type AbstractEvolutionState end
abstract type AbstractNonlinearSolver end
abstract type AbstractTimestepController end
abstract type AbstractTerminationCondition end
abstract type AbstractRetryPolicy end
abstract type AbstractRemeshPolicy end


Base.@kwdef struct StellarEvolver{TSOLVER<:AbstractNonlinearSolver, TTIMESTEP<:AbstractTimestepController, TTERMINATION<:AbstractTerminationCondition, 
    TRETRY<:AbstractRetryPolicy,TREMESH<:AbstractRemeshPolicy,} <: AbstractEvolver
    solver::TSOLVER
    timestep::TTIMESTEP
    termination::TTERMINATION
    retry::TRETRY
    remesh::TREMESH
end

Base.@kwdef struct OneZoneEvolver{TSOLVER<:AbstractNonlinearSolver, TTIMESTEP<:AbstractTimestepController, TTERMINATION<:AbstractTerminationCondition, 
    TRETRY<:AbstractRetryPolicy,} <: AbstractEvolver
    solver::TSOLVER
    timestep::TTIMESTEP
    termination::TTERMINATION
    retry::TRETRY
end

Base.@kwdef mutable struct StellarEvolutionState <: AbstractEvolutionState
    timestep_retry_count::Int = 0
    base_timestep::Float64 = 0.0
    use_remeshing::Bool = true
    starting_new_step::Bool = true
end

Base.@kwdef mutable struct OneZoneEvolutionState <: AbstractEvolutionState
    timestep_retry_count::Int = 0
    base_timestep::Float64 = 0.0
    starting_new_step::Bool = true
end

struct StepResult
    converged::Bool 
    reason::Symbol
    iterations::Int
end 



#intended user use 
# stellar_evolver = StellarEvolver(
#     solver,
#     timestep_controller,
#     termination_conditions,
#     retry_policy,
#     remesh_policy,
# )

# evolve!(sm, stellar_evolver)