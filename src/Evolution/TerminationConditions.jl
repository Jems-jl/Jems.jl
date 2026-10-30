Base.@kwdef struct StellarTermination <: AbstractTerminationCondition
    maximum_model_number::Union{Nothing,Int} = nothing
    maximum_age::Union{Nothing,Float64} = nothing
    maximum_center_temperature::Union{Nothing,Float64} = nothing
    minimum_center_hydrogen::Union{Nothing,Float64} = nothing
    extra_conditions::Tuple = ()

end 

function termination_reason(condition::StellarTermination, sm::StellarModel, evolver, state)
    if condition.maximum_model_number !== nothing && sm.props.model_number > condition.maximum_model_number
        return :maximum_model_number
    end

    if condition.maximum_center_temperature !== nothing && exp(get_value(sm.props.lnT[1])) > condition.maximum_center_temperature
        return :maximum_center_temperature
    end

    center_hydrogen = get_value(sm.props.xa[1, sm.network.xa_index[:H1]])

    if condition.minimum_center_hydrogen !== nothing && center_hydrogen < condition.minimum_center_hydrogen
        return :minimum_center_hydrogen
    end
    
    if condition.maximum_age !== nothing && sm.props.time > condition.maximum_age
        return :maximum_time 
    end
    for extra_condition in condition.extra_conditions
        reason = termination_reason( extra_condition, sm, evolver, state,)
        reason === nothing || return reason
    end

    return nothing
end 