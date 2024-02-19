# """
#     set_end_step_info(sm::StellarModel)

# Sets the StellarStepInfo `si`` from current state of the StellarModel `sm`.
# """
# function set_step_info!(sm::StellarModel, si::StellarModels.StellarStepInfo)
#     si.model_number = sm.model_number
#     si.time = sm.time
#     si.dt = sm.dt

#     si.nz = sm.props.nz
#     si.mstar = sm.mstar

#     Threads.@threads for i = 1:(sm.props.nz)
#         si.m[i] = sm.m[i]
#         si.dm[i] = sm.dm[i]

#         si.lnT[i] = sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnT]]
#         si.L[i] = sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lum]]
#         si.lnρ[i] = sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnρ]]
#         si.lnr[i] = sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnr]]

#         xa = view(sm.ind_vars, (i * sm.nvars - sm.network.nspecies + 1):(i * sm.nvars))
#         si.X[i] = xa[sm.network.xa_index[:H1]]  # can later include H2 as well.
#         si.Y[i] = xa[sm.network.xa_index[:He4]]  # can later include He3 as well.

#         set_EOS_resultsTρ!(sm.eos, si.eos_res[i], si.lnT[i], si.lnρ[i], xa, sm.network.species_names)

#         si.lnP[i] = log(si.eos_res[i].P)
#         for k = 1:sm.nvars
#             si.ind_vars[(i - 1) * sm.nvars + k] = sm.ind_vars[(i - 1) * sm.nvars + k]
#         end
#     end
# end

"""
    cycle_props!(sm::StellarModel)

Moves the model properties of the StellarModel `sm` over one state:
start_step_props -> props -> prv_step_props -> start_step_props
"""
function cycle_props!(sm::StellarModel)
    temp_props = sm.prv_step_props
    sm.prv_step_props = sm.props
    sm.props = sm.start_step_props
    sm.start_step_props = temp_props
end

"""
    uncycle_props!(sm::StellarModel)

Moves the model properties of the StellarModel `sm` back one state:
start_step_props <- props <- prv_step_props <- start_step_props
"""
function uncycle_props!(sm::StellarModel)
    temp_props = sm.props
    sm.props = sm.prv_step_props
    sm.prv_step_props = sm.start_step_props
    sm.start_step_props = temp_props
end

"""
    get_dt_next(sm::StellarModel)

Computes the timestep of the next evolutionary step to be taken by the StellarModel `sm` by considering all timestep
controls (`sm.opt.timestep`).
"""
function get_dt_next(sm::StellarModel)
    dt_next = sm.props.dt  # this it calculated at end of step, so props.dt is the dt we used to do this step
        
    Rsurf = exp(get_cell_value(sm.props.lnr[sm.props.nz]))
    Rsurf_old = exp(get_cell_value(sm.prv_step_props.lnr[sm.prv_step_props.nz]))
    ΔR_div_R = abs(Rsurf - Rsurf_old) / Rsurf

    Tc = exp(get_cell_value(sm.props.lnT[sm.props.nz]))
    Tc_old = exp(get_cell_value(sm.prv_step_props.lnT[sm.prv_step_props.nz]))
    ΔTc_div_Tc = abs(Tc - Tc_old) / Tc

    X = get_cell_value(sm.props.xa[sm.props.nz, sm.network.xa_index[:H1]])
    Xold = get_cell_value(sm.prv_step_props.xa[sm.prv_step_props.nz, sm.network.xa_index[:H1]])
    ΔX = abs(X - Xold) / (X)

    dt_nextR = dt_next * sm.opt.timestep.delta_R_limit / ΔR_div_R
    dt_nextTc = dt_next * sm.opt.timestep.delta_Tc_limit / ΔTc_div_Tc
    dt_nextX = dt_next * sm.opt.timestep.delta_Xc_limit / ΔX

    min_dt = dt_next * sm.opt.timestep.dt_max_decrease
    dt_next = min(sm.opt.timestep.dt_max_increase * dt_next, dt_nextR, dt_nextTc, dt_nextX)
    dt_next = max(dt_next, min_dt)
    return dt_next
end

"""
    do_evolution_loop(sm::StellarModel)

Performs the main evolutionary loop of the input StellarModel `sm`. It continues taking steps until one of the
termination criteria is reached (defined in `sm.opt.termination`).
"""
function do_evolution_loop!(sm::StellarModel)
    # before loop actions
    StellarModels.create_output_files!(sm)
    StellarModels.update_stellar_model_properties!(sm, sm.props)  # set the initial condition as the result of a previous phantom step
    dt_factor = 1.0  # this is changed during retries to lower the timestep
    retry_count = 0

    # evolution loop, be sure to have sensible termination conditions or this will go on forever!
    while true
        cycle_props!(sm)  # move props of previous step to prv_step_props of current step

        # remeshing
        if sm.opt.remesh.do_remesh
            sm = StellarModels.remesher!(sm)
        end

        # MF: this is unnecessary atm I think, we never refer to start_step_props
        StellarModels.update_stellar_model_properties!(sm, sm.start_step_props)  # save start_step_props before we attempt any newton solver

        sm.solver_data.newton_iters = 0
        max_steps = sm.opt.solver.newton_max_iter
        if (sm.model_number == 0)
            max_steps = sm.opt.solver.newton_max_iter_first_step
        end

        exit_evolution = false
        retry_step = false
        # step loop
        for i = 1:max_steps
            StellarModels.update_stellar_model_properties!(sm, sm.props)

            eval_jacobian_eqs!(sm)  # heavy lifting happens here!
            thomas_algorithm!(sm)  # here as well

            corr = @view sm.solver_data.solver_corr[1:sm.nvars*sm.props.nz]
            real_max_corr = maximum(corr)

            # scale surface correction to prevent negative surface luminosity
            # if correction will produce negative L, scale it so L is halved
            corr_lum_surf = corr[sm.nvars * (sm.props.nz - 1) + sm.vari[:lum]]
            lum_surf = sm.ind_vars[sm.nvars * (sm.props.nz - 1) + sm.vari[:lum]]
            if lum_surf + corr_lum_surf < 0.0
                corr *= (-0.1 * lum_surf / corr_lum_surf)
            end

            # scale correction
            if sm.model_number == 0
                corr .*= min(1, sm.opt.solver.initial_model_scale_max_correction / maximum(corr))
            else
                corr .*= min(1, sm.opt.solver.scale_max_correction / maximum(corr))
            end
            # if i % 50 == 0
            #     @show i, maximum(corr), real_max_corr, maximum(sm.solver_data.eqs_numbers)
            # end

            # applying correction
            sm.ind_vars[1:sm.nvars*sm.props.nz] .+= corr[1:sm.nvars*sm.props.nz]
            if real_max_corr < 1e-10
                if sm.model_number == 0
                    println("Found first model")
                end
                break  # successful, break the step loop
            end
            if i == max_steps
                if retry_count > 10
                    exit_evolution = true
                    println("Too many retries, ending simulation")
                else
                    retry_count = retry_count + 1
                    retry_step = true
                    println("Failed to converge step $(sm.model_number) with timestep $(sm.props.dt/SECYEAR), retrying")
                end
            end
            sm.solver_data.newton_iters = i
        end

        if retry_step
            dt_factor *= dt_retry_decrease
            # adapt dt for coming step
            sm.dt *= dt_factor
            uncycle_props!(sm)  # reset props to what prv_step_props contains
            continue
        else
            dt_factor = 1.0
        end

        if (exit_evolution)
            println("Terminating evolution")
            break
        end

        # step must be successful at this point
        retry_count = 0

        # increment age and model number since we accept the step.
        sm.time += sm.props.dt
        sm.model_number += 1

        # write state in sm.props and potential history/profiles.
        StellarModels.update_stellar_model_properties!(sm, sm.props)
        # StellarModels.write_data(sm)
        StellarModels.write_terminal_info(sm)

        if sm.opt.plotting.do_plotting && sm.model_number == 1
            Plotting.init_plots!(sm)
        elseif sm.opt.plotting.do_plotting && sm.model_number % sm.opt.plotting.plotting_interval == 0
            Plotting.update_plotting!(sm)
        end

        #@show sm.model_number, sm.esi.lnP[1], sm.esi.lnP[2], sm.esi.lnP[sm.props.nz-1], sm.esi.lnP[sm.props.nz]
        #@show sm.model_number, sm.esi.lnT[1], sm.esi.lnT[2], sm.esi.lnT[sm.props.nz-1], sm.esi.lnT[sm.props.nz]
        #@show sm.dm[1], sm.dm[2], sm.dm[3]
        #@show sum(sm.dm[1:sm.props.nz])

        # check termination conditions
        if (sm.model_number > sm.opt.termination.max_model_number)
            println("Reached maximum model number")
            break
        end
        if (exp(get_cell_value(sm.props.lnT[1])) > sm.opt.termination.max_center_T)
            println("Reached maximum central temperature")
            break
        end

        # get dt for coming step
        sm.dt = get_dt_next(sm)
    end
    if sm.opt.plotting.do_plotting
        Plotting.end_of_evolution(sm)
    end
    StellarModels.close_output_files!(sm)
    return sm
end
