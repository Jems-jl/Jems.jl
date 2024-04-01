
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
    do_evolution_loop(sm::StellarModel)

Performs the main evolutionary loop of the input StellarModel `sm`. It continues taking steps until one of the
termination criteria is reached (defined in `sm.opt.termination`).
"""
function do_evolution_loop!(sm::StellarModel)
    # before loop actions
    StellarModels.create_output_files!(sm)
    StellarModels.evaluate_stellar_model_properties!(sm, sm.props)  # set the initial condition as the result of a previous phantom step
    retry_count = 0

    # evolution loop, be sure to have sensible termination conditions or this will go on forever!
    while true
        cycle_props!(sm)  # move props of previous step to prv_step_props of current step

        # either remesh, or copy over from prv_step_props
        if sm.opt.remesh.do_remesh
            StellarModels.remesher!(sm)
        else
            StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)
        end
        # time derivatives in the equations use the remeshed info, so
        # save start_step_props before we attempt any newton solver
        StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
        StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
        sm.start_step_props.dt = sm.prv_step_props.dt_next  # dt of this step becomes dt_next of previous

        # step loop
        sm.solver_data.newton_iters = 0
        max_steps = sm.opt.solver.newton_max_iter
        if (sm.start_step_props.model_number == 0)
            max_steps = sm.opt.solver.newton_max_iter_first_step
        end

        exit_evolution = false
        retry_step = false
        StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
        StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

        corr = @view sm.solver_data.solver_corr[1:sm.nvars*sm.props.nz]
        equs = @view sm.solver_data.eqs_numbers[1:sm.nvars*sm.props.nz]

        # evaluate the equations for the first step
        StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
        eval_jacobian_eqs!(sm)  # heavy lifting happens here!
        for i = 1:max_steps
            thomas_algorithm!(sm)  # here as well

            (abs_max_corr, i_corr) = findmax(abs, corr)
            signed_max_corr = corr[i_corr]
            corr_nz = i_corr÷sm.nvars + 1
            corr_equ = i_corr%sm.nvars
            rel_corr = abs_max_corr/eps(sm.props.ind_vars[i_corr])

            # scale correction
            if sm.props.model_number == 0
                correction_multiplier = min(1.0, sm.opt.solver.initial_model_scale_max_correction / abs_max_corr)
            else
                correction_multiplier = min(1.0, sm.opt.solver.scale_max_correction / abs_max_corr)
            end
            if correction_multiplier < 1
                corr .*= correction_multiplier
            end

            # first try applying correction and see if it would give negative luminosity
            sm.props.ind_vars[1:sm.nvars*sm.props.nz] .+= corr[1:sm.nvars*sm.props.nz]
            sm.solver_data.newton_iters = i

            # evaluate the equations after correction and get residuals
            try
                StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
                eval_jacobian_eqs!(sm)  # heavy lifting happens here!

                (max_res, i_res) = findmax(abs, equs)
                res_nz = i_res÷sm.nvars + 1
                res_equ = i_res%sm.nvars

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
                    break  # successful, break the step loop
                end
            catch e
                if isa(e, InterruptException)
                    throw(e)
                end
                println("Error while evaluating equations")
                showerror(stdout, e)
                retry_step = true
            end
            #if not, determine if we give up or retry
            if i == max_steps
                if retry_count > 10
                    exit_evolution = true
                    println("Too many retries, ending simulation")
                else
                    retry_count = retry_count + 1
                    retry_step = true
                    println("Failed to converge step $(sm.props.model_number) with timestep $(sm.props.dt/SECYEAR), retrying")
                end
            end
        end

        if retry_step
            uncycle_props!(sm)  # reset props to what prv_step_props contains, ie mimic state at end of previous step
            sm.props.dt_next *= sm.opt.timestep.dt_retry_decrease # adapt dt
            continue  # go back to top of evolution loop
        end

        if (exit_evolution)
            println("Terminating evolution")
            break
        end

        # step must be successful at this point
        retry_count = 0

        # increment age and model number since we accept the step.
        sm.props.time += sm.props.dt
        sm.props.model_number += 1

        # write state in sm.props and potential history/profiles.
        StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
        StellarModels.write_data(sm)
        StellarModels.write_terminal_info(sm)

        if sm.opt.plotting.do_plotting && sm.props.model_number == 1
            Plotting.init_plots!(sm)
        elseif sm.opt.plotting.do_plotting && sm.props.model_number % sm.opt.plotting.plotting_interval == 0
            Plotting.update_plotting!(sm)
        end

        # check termination conditions
        if (sm.props.model_number > sm.opt.termination.max_model_number)
            StellarModels.write_terminal_info(sm; now=true)
            println("Reached maximum model number")
            break
        end
        if (exp(get_value(sm.props.lnT[1])) > sm.opt.termination.max_center_T)
            StellarModels.write_terminal_info(sm; now=true)
            println("Reached maximum central temperature")
            break
        end

        # get dt for coming step
        sm.props.dt_next = get_dt_next(sm)
    end
    if sm.opt.plotting.do_plotting
        Plotting.end_of_evolution(sm)
    end
    StellarModels.close_output_files!(sm)
end
