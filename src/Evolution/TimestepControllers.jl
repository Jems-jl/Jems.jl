struct StellarTimestepController <: AbstractTimestepController
end
struct OneZoneTimestepController <: AbstractTimestepController
end
struct UserTimestepController <: AbstractTimestepController
end 

function get_dt_next(timestep::StellarTimestepController, sm::StellarModel, evolver, state)
    dt_next = sm.props.dt  # this it calculated at end of step, so props.dt is the dt we used to do this step
        
    Rsurf = exp(get_value(sm.props.lnr[sm.props.nz]))
    Rsurf_old = exp(get_value(sm.prv_step_props.lnr[sm.prv_step_props.nz]))
    ΔR_div_R = abs(Rsurf - Rsurf_old) / Rsurf

    Tsurf = exp(get_value(sm.props.lnT[sm.props.nz]))
    Tsurf_old = exp(get_value(sm.prv_step_props.lnT[sm.prv_step_props.nz]))
    ΔTc_div_Tc = abs(Tsurf - Tsurf_old) / Tsurf

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

get_next_timestep(timestep::OneZoneTimestepController,::OneZone,evolver,state)
