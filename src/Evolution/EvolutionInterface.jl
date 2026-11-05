abstract type AbstractEvolver end
abstract type AbstractNonlinearSolver end

struct StellarEvolver{TSOLVER<:AbstractNonlinearSolver} <: AbstractEvolver
    solver::TSOLVER
    max_timestep_retries::Int
    extra_termination_conditions::Tuple
end

struct OneZoneEvolver{TSOLVER<:AbstractNonlinearSolver} <: AbstractEvolver
    solver::TSOLVER
    max_timestep_retries::Int
    extra_termination_conditions::Tuple
end

Base.@kwdef mutable struct StellarEvolutionState
    timestep_retry_count::Int = 0
    base_timestep::Float64 = 0.0
    use_remeshing::Bool = true
    starting_new_step::Bool = true
end

Base.@kwdef mutable struct OneZoneEvolutionState
    timestep_retry_count::Int = 0
    base_timestep::Float64 = 0.0
    starting_new_step::Bool = true
end

struct StepResult
    converged::Bool
    reason::Symbol
    iterations::Int
end
