function get_dt_next(oz::OneZone)
    dt_next = oz.props.dt  # this it calculated at end of step, so props.dt is the dt we used to do this step

    X = get_value(oz.props.xa[oz.network.xa_index[:H1]])
    Xold = get_value(oz.prv_step_props.xa[oz.network.xa_index[:H1]])
    ΔX = abs(X - Xold)

    dt_nextX = dt_next * oz.opt.timestep.delta_Xc_limit / ΔX

    min_dt = dt_next * oz.opt.timestep.dt_max_decrease
    dt_next = min(oz.opt.timestep.dt_max_increase * dt_next, dt_nextX)
    dt_next = max(dt_next, min_dt)
    dt_next = min(dt_next, oz.opt.timestep.max_dt * SECYEAR)
    return dt_next
end

function do_one_zone_burn!(oz::OneZone)
    # before loop actions
    StellarModels.create_output_files!(oz)
    StellarModels.evaluate_model_properties!(oz, oz.props)  # set the initial condition as the result of a previous phantom step
    retry_count = 0

    # evolution loop, be sure to have sensible termination conditions or this will go on forever!
    while true
        StellarModels.cycle_props!(oz)  # move props of previous step to prv_step_props of current step

        # step loop
        oz.solver_data.newton_iters = 0
        max_steps = oz.opt.solver.newton_max_iter
        if (oz.props.model_number == 0)
            max_steps = oz.opt.solver.newton_max_iter_first_step
        end

        exit_evolution = false
        retry_step = false
        oz.props = deepcopy(oz.prv_step_props)
        oz.props.dt = oz.prv_step_props.dt_next  # dt of this step becomes dt_next of previous

        corr = oz.solver_data.solver_corr
        equs = oz.solver_data.eqs_numbers

        # evaluate the equations for the first step
        eval_jacobian_eqs!(oz)  # heavy lifting happens here!
        for i = 1:max_steps
            thomas_algorithm!(oz)  # here as well

            (abs_max_corr, i_corr) = findmax(abs, corr)
            signed_max_corr = corr[i_corr]
            corr_nz = i_corr ÷ oz.nvars + 1
            corr_equ = i_corr % oz.nvars
            rel_corr = abs_max_corr / eps(oz.props.ind_vars[i_corr])

            # scale correction
            if oz.props.model_number == 0
                correction_multiplier = min(1.0, oz.opt.solver.initial_model_scale_max_correction / abs_max_corr)
            else
                correction_multiplier = min(1.0, oz.opt.solver.scale_max_correction / abs_max_corr)
            end
            if correction_multiplier < 1
                corr .*= correction_multiplier
            end

            # first try applying correction and see if it would give negative luminosity
            oz.props.ind_vars .+= corr
            oz.solver_data.newton_iters = i

            # evaluate the equations after correction and get residuals
            try
                StellarModels.evaluate_model_properties!(oz, oz.props)
                eval_jacobian_eqs!(oz)  # heavy lifting happens here!

                (max_res, i_res) = findmax(abs, equs)
                res_nz = i_res ÷ oz.nvars + 1
                res_equ = i_res % oz.nvars

                #reporting
                if oz.opt.solver.report_solver_progress &&
                   i % oz.opt.solver.solver_progress_iter == 0
                    @show oz.props.model_number, i, rel_corr, signed_max_corr, corr_nz, corr_equ, max_res, res_nz,
                          res_equ
                end
                #check if tolerances are satisfied
                if rel_corr < oz.opt.solver.relative_correction_tolerance &&
                   max_res < oz.opt.solver.maximum_residual_tolerance
                    if oz.props.model_number == 0
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
            # if not, determine if we give up or retry
            if i == max_steps
                if retry_count > 10
                    exit_evolution = true
                    println("Too many retries, ending simulation")
                else
                    retry_count = retry_count + 1
                    retry_step = true
                    println("Failed to converge step $(oz.props.model_number) with timestep $(oz.props.dt/SECYEAR), retrying")
                end
            end
        end

        if retry_step
            StellarModels.uncycle_props!(oz)  # reset props to what prv_step_props contains, ie mimic state at end of previous step
            oz.props.dt_next *= oz.opt.timestep.dt_retry_decrease # adapt dt
            continue  # go back to top of evolution loop
        end

        if (exit_evolution)
            println("Terminating evolution")
            break
        end

        # step must be successful at this point
        retry_count = 0

        # increment age and model number since we accept the step.
        oz.props.time += oz.props.dt
        oz.props.model_number += 1

        # write state in oz.props and potential history/profiles.
        StellarModels.evaluate_model_properties!(oz, oz.props)
        StellarModels.write_data(oz)
        StellarModels.write_terminal_info(oz)

        if oz.opt.plotting.do_plotting && oz.props.model_number == 1
            Plotting.init_plots!(oz)
        elseif oz.opt.plotting.do_plotting && oz.props.model_number % oz.opt.plotting.plotting_interval == 0
            Plotting.update_plotting!(oz)
        end

        # check termination conditions
        if (oz.props.model_number > oz.opt.termination.max_model_number)
            # StellarModels.write_terminal_info(oz; now=true)
            println("Reached maximum model number")
            break
        end

        # get dt for coming step
        oz.props.dt_next = get_dt_next(oz)
    end
    if oz.opt.plotting.do_plotting
        Plotting.end_of_evolution(oz)
    end
    StellarModels.shut_down_IO!(oz)
end
