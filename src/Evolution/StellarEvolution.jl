
# """
#     get_dt_next(sm::StellarModel)

# Computes the timestep of the next evolutionary step to be taken by the StellarModel `sm` by considering all timestep
# controls (`sm.opt.timestep`).
# """
# function get_dt_next(sm::StellarModel)
#     dt_next = sm.props.dt  # this it calculated at end of step, so props.dt is the dt we used to do this step
        
#     Rsurf = exp(get_value(sm.props.lnr[sm.props.nz]))
#     Rsurf_old = exp(get_value(sm.prv_step_props.lnr[sm.prv_step_props.nz]))
#     ΔR_div_R = abs(Rsurf - Rsurf_old) / Rsurf

#     Tc = exp(get_value(sm.props.lnT[sm.props.nz]))
#     Tc_old = exp(get_value(sm.prv_step_props.lnT[sm.prv_step_props.nz]))
#     ΔTc_div_Tc = abs(Tc - Tc_old) / Tc

#     X = get_value(sm.props.xa[1, sm.network.xa_index[:H1]])
#     Xold = get_value(sm.prv_step_props.xa[1, sm.network.xa_index[:H1]])
#     ΔX = abs(X - Xold)

#     dt_nextR = dt_next * sm.opt.timestep.delta_R_limit / ΔR_div_R
#     dt_nextTc = dt_next * sm.opt.timestep.delta_Tc_limit / ΔTc_div_Tc
#     dt_nextX = dt_next * sm.opt.timestep.delta_Xc_limit / ΔX

#     min_dt = dt_next * sm.opt.timestep.dt_max_decrease
#     dt_next = min(sm.opt.timestep.dt_max_increase * dt_next, dt_nextR, dt_nextTc, dt_nextX)
#     dt_next = max(dt_next, min_dt)
#     return dt_next
# end

# """
#     compute_starting_model_properties(sm::StellarModel)

# Computes all stellar model properties based on the current state of the star.
# This sets the model up to start the evolution loop, or test the evaluation of
# the equations and the linear solver.
# """
# function compute_starting_model_properties!(sm::StellarModel)
#     StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#     StellarModels.cycle_props!(sm);
#     StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
#     StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)
#     StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
#     StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
#     StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)
# end

# """
#     do_evolution_loop(sm::StellarModel)

# Performs the main evolutionary loop of the input StellarModel `sm`. It continues taking steps until one of the
# termination criteria is reached (defined in `sm.opt.termination`).
# """
# function do_evolution_loop!(sm::StellarModel; plotter::TPLOTTER = Plotting.NullPlotter(),
#                                 max_retries_in_a_row = 10) where{TPLOTTER<:AbstractPlotter}
#     # before loop actions
#     TNUMBER = eltype(sm.props.ind_vars)  # determine type of numbers
#     StellarModels.create_output_files!(sm, TNUMBER)
#     compute_starting_model_properties!(sm)
#     retry_count = 0

#     # evolution loop, be sure to have sensible termination conditions or this will go on forever!
#     while true
#         StellarModels.cycle_props!(sm)  # move props of previous step to prv_step_props of current step

#         # either remesh, or copy over from prv_step_props
#         if sm.opt.remesh.do_remesh
#             StellarModels.remesher!(sm)
#         else
#             StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)
#         end
#         # time derivatives in the equations use the remeshed info, so
#         # save start_step_props before we attempt any newton solver
#         StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
#         StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
#         sm.start_step_props.dt = sm.prv_step_props.dt_next  # dt of this step becomes dt_next of previous

#         # step loop
#         sm.solver_data.newton_iters = 0
#         max_steps = sm.opt.solver.newton_max_iter
#         if (sm.start_step_props.model_number == 0)
#             max_steps = sm.opt.solver.newton_max_iter_first_step
#         end

#         exit_evolution = false
#         retry_step = false
#         StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
#         StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

#         corr = @view sm.solver_data.solver_corr[1:sm.nvars*sm.props.nz]
#         equs = @view sm.solver_data.eqs_numbers[1:sm.nvars*sm.props.nz]

#         # evaluate the equations for the first step
#         StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#         eval_jacobian_eqs!(sm)  # heavy lifting happens here!
#         for i = 1:max_steps
#             block_tridiagonal_solver!(sm, sm.solver_data)  # here as well

#             # (abs_max_corr, i_corr) = findmax(abs, corr)
#             # signed_max_corr = corr[i_corr]
#             # corr_nz = i_corr÷sm.nvars + 1
#             # corr_equ = i_corr%sm.nvars
#             # rel_corr = abs_max_corr/eps(sm.props.ind_vars[i_corr])

#             # scale correction
#             if sm.props.model_number == 0
#                 correction_limit = sm.opt.solver.initial_model_scale_max_correction
#             else
#                 correction_limit = sm.opt.solver.scale_max_correction
#             end
#             correction_multiplier = 1.0
#             for j in 1:sm.nvars
#                 sub_corr = @view sm.solver_data.solver_corr[j:sm.nvars:(sm.nvars*(sm.props.nz-1)+j)]
#                 max_sub_corr = maximum(abs, sub_corr)
#                 if sm.var_scaling[j] == :log || sm.var_scaling[j] == :unity
#                     correction_multiplier = min(correction_multiplier,
#                                                     correction_limit/max_sub_corr)
#                 elseif sm.var_scaling[j] == :maxval
#                     sub_ind_vars = @view sm.props.ind_vars[j:sm.nvars:(sm.nvars*(sm.props.nz-1)+j)]
#                     max_var_value = maximum(abs, sub_ind_vars)
#                     correction_multiplier = min(correction_multiplier,
#                                                     max_var_value*correction_limit/max_sub_corr)
#                 end
#             end
            
#             if correction_multiplier < 1.0
#                 corr .*= correction_multiplier
#             end

#             gamma_j = sm.vari[:gamma_turb]
#             gamma_min = 0.0
#             boundary_fraction = 0.9

#             for k in 1:sm.props.nz
#                 idx = (k - 1) * sm.nvars + gamma_j

#                 gamma_old = sm.props.ind_vars[idx]
#                 delta_gamma = corr[idx]

#                 if gamma_old <= gamma_min
#                     # This should only be needed when restarting from an old profile
#                     # that already contains nonpositive gamma.
#                     gamma_old = gamma_min
#                     sm.props.ind_vars[idx] = gamma_min
#                 end

#                 # Limit only a negative correction that would cross gamma_min.
#                 # Stay 10% of the old distance away from the boundary rather than
#                 # landing exactly on gamma_min.
#                 maximum_downward_step =
#                     boundary_fraction * (gamma_old - gamma_min)

#                 corr[idx] = max(delta_gamma, -maximum_downward_step)
#             end
#             (abs_max_corr, i_corr) = findmax(abs, corr)
#             signed_max_corr = corr[i_corr]

#             corr_nz = div(i_corr - 1, sm.nvars) + 1
#             corr_equ = mod(i_corr - 1, sm.nvars) + 1

#             # Retain the original convergence definition.
#             rel_corr = abs_max_corr / eps(sm.props.ind_vars[i_corr])
#             sm.props.ind_vars[1:sm.nvars*sm.props.nz] .+= corr[1:sm.nvars*sm.props.nz]
#             sm.solver_data.newton_iters = i


#             # evaluate the equations after correction and get residuals
#             try
#                 StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#                 eval_jacobian_eqs!(sm)  # heavy lifting happens here!

#                 (max_res, i_res) = findmax(abs, equs)
#                 # res_nz = i_res÷sm.nvars + 1
#                 # res_equ = i_res%sm.nvars
#                 res_nz = div(i_res - 1, sm.nvars) + 1
#                 res_equ = mod(i_res - 1, sm.nvars) + 1

#                 #reporting
#                 if sm.opt.solver.report_solver_progress &&
#                     i % sm.opt.solver.solver_progress_iter == 0
#                     @show sm.props.model_number, i, rel_corr, signed_max_corr, corr_nz, corr_equ, max_res, res_nz, res_equ
#                 end
#                 #check if tolerances are satisfied
#                 if rel_corr < sm.opt.solver.relative_correction_tolerance &&
#                         max_res < sm.opt.solver.maximum_residual_tolerance
#                     if sm.props.model_number == 0
#                         println("Found first model")
#                     end
#                     break  # successful, break the step loop
#                 end
#             catch e
#                 if isa(e, InterruptException)
#                     throw(e)
#                 end
#                 println("Error while evaluating equations")
#                 showerror(stdout, e)
#                 retry_step = true
#             end
#             #if not, determine if we give up or retry
#             if i == max_steps
#                 if retry_count > max_retries_in_a_row
#                     exit_evolution = true
#                     println("Too many retries, ending simulation")
#                 else
#                     retry_count = retry_count + 1
#                     retry_step = true
#                     if sm.opt.solver.report_retries
#                         println("Failed to converge step $(sm.props.model_number) with timestep $(sm.props.dt/SECYEAR), retrying")
#                     end
#                 end
#             end
#         end

        # n = sm.nvars * sm.props.nz
        # x = @view sm.props.ind_vars[1:n]
        # x0 = similar(x)

        # sufficient_decrease = 1e-8
        # max_backtracks = 16

        # for i = 1:max_steps
        #     # The current residual and Jacobian are already valid.
        #     # Save the baseline BEFORE the linear solver modifies its arrays.
        #     copyto!(x0, x)
        #     R0 = maximum(abs, equs)
        #     isfinite(R0) || error("Nonfinite baseline residual")

        #     block_tridiagonal_solver!(sm, sm.solver_data)

        #     # Keep convergence diagnostics based on the raw Newton correction.
        #     abs_max_corr, i_corr = findmax(abs, corr)
        #     signed_max_corr = corr[i_corr]
        #     corr_nz = div(i_corr - 1, sm.nvars) + 1
        #     corr_equ = mod(i_corr - 1, sm.nvars) + 1
        #     rel_corr = abs_max_corr / eps(x0[i_corr])

        #     # Existing correction limits.
        #     correction_limit = if sm.props.model_number == 0
        #         sm.opt.solver.initial_model_scale_max_correction
        #     else
        #         sm.opt.solver.scale_max_correction
        #     end

        #     correction_multiplier = 1.0

        #     for j in 1:sm.nvars
        #         sub_corr = @view corr[j:sm.nvars:n]
        #         max_sub_corr = maximum(abs, sub_corr)

        #         # Avoid division by zero for an unchanged variable.
        #         max_sub_corr == 0 && continue

        #         if sm.var_scaling[j] == :log ||
        #         sm.var_scaling[j] == :unity
        #             correction_multiplier = min(
        #                 correction_multiplier,
        #                 correction_limit / max_sub_corr,
        #             )
        #         elseif sm.var_scaling[j] == :maxval
        #             sub_vars = @view x0[j:sm.nvars:n]
        #             max_var_value = maximum(abs, sub_vars)

        #             correction_multiplier = min(
        #                 correction_multiplier,
        #                 max_var_value * correction_limit / max_sub_corr,
        #             )
        #         end
        #     end

        #     # p remains fixed throughout this line search.
        #     p = correction_multiplier .* corr

        #     accepted = false
        #     converged = false
        #     accepted_alpha = 0.0
        #     max_res = R0

        #     try
        #         valid_direction = all(isfinite, p) &&
        #             isfinite(correction_multiplier) &&
        #             0 < correction_multiplier <= 1

        #         if valid_direction
        #             for attempt in 0:max_backtracks
        #                 alpha = 0.5^attempt

        #                 # attempt=0 is the usual full limited Newton step.
        #                 # Smaller fractions are tried only after rejection.
        #                 x .= x0 .+ alpha .* p

        #                 Rtrial = try
        #                     StellarModels.evaluate_stellar_model_properties!(
        #                         sm, sm.props,
        #                     )
        #                     eval_jacobian_eqs!(sm)
        #                     maximum(abs, equs)
        #                 catch err
        #                     if err isa DomainError || err isa OverflowError
        #                         Inf
        #                     else
        #                         rethrow()
        #                     end
        #                 end

        #                 trial_converged = isfinite(Rtrial) &&
        #                     rel_corr <
        #                         sm.opt.solver.relative_correction_tolerance &&
        #                     Rtrial <
        #                         sm.opt.solver.maximum_residual_tolerance

        #                 target = (
        #                     1 - sufficient_decrease *
        #                         alpha * correction_multiplier
        #                 ) * R0

        #                 sufficient_improvement = isfinite(Rtrial) &&
        #                     Rtrial < R0 &&
        #                     Rtrial <= target

        #                 if trial_converged || sufficient_improvement
        #                     accepted = true
        #                     converged = trial_converged
        #                     accepted_alpha = alpha
        #                     max_res = Rtrial
        #                     break
        #                 end
        #             end
        #         end
        #         finally
        #             if !accepted
        #                 # Restore a consistent baseline on failure or exception.
        #                 copyto!(x, x0)
        #                 StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
        #                 eval_jacobian_eqs!(sm)
        #             end
        #         end

        #         if !accepted
        #             println(
        #                 "Backtracking failed: model=", sm.props.model_number,
        #                 ", iteration=", i,
        #                 ", residual=", R0,
        #             )

        #             if retry_count >= max_retries_in_a_row
        #                 retry_step = false
        #                 exit_evolution = true
        #                 println("Too many retries, ending simulation")
        #             else
        #                 retry_count += 1
        #                 retry_step = true
        #             end

        #             break
        #         end

        #         # x, properties, residuals, and Jacobian ALREADY describe
        #         # the accepted state. Do not apply or evaluate it again.
        #         corr .= accepted_alpha .* p
        #         sm.solver_data.newton_iters = i

        #         _, i_res = findmax(abs, equs)
        #         res_nz = div(i_res - 1, sm.nvars) + 1
        #         res_equ = mod(i_res - 1, sm.nvars) + 1

        #         if accepted_alpha < 1.0
        #             println(
        #                 "Backtracking accepted: model=", sm.props.model_number,
        #                 ", iteration=", i,
        #                 ", alpha=", accepted_alpha,
        #                 ", residual=", max_res,
        #             )
        #         end

        #         if sm.opt.solver.report_solver_progress &&
        #         i % sm.opt.solver.solver_progress_iter == 0
        #             @show sm.props.model_number, i, rel_corr, signed_max_corr,
        #                 corr_nz, corr_equ, max_res, res_nz, res_equ
        #         end

        #         if converged
        #             if sm.props.model_number == 0
        #                 println("Found first model")
        #             end
        #             break
        #         end

        #         if i == max_steps
        #             if retry_count >= max_retries_in_a_row
        #                 retry_step = false
        #                 exit_evolution = true
        #                 println("Too many retries, ending simulation")
        #             else
        #                 retry_count += 1
        #                 retry_step = true

        #                 if sm.opt.solver.report_retries
        #                     println(
        #                         "Newton iteration limit reached: model=",
        #                         sm.props.model_number,
        #                         ", dt=", sm.props.dt / SECYEAR, " years",
        #                     )
        #                 end
        #             end
        #         end
        #     end

#         if retry_step
#             StellarModels.uncycle_props!(sm)  # reset props to what prv_step_props contains, ie mimic state at end of previous step
#             sm.props.dt_next *= sm.opt.timestep.dt_retry_decrease # adapt dt
#             continue  # go back to top of evolution loop
#         end

#         if (exit_evolution)
#             println("Terminating evolution")
#             break
#         end

#         # step must be successful at this point
#         retry_count = 0

#         # increment age and model number since we accept the step.
#         sm.props.time += sm.props.dt
#         sm.props.model_number += 1

#         # write state in sm.props and potential history/profiles.
#         StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#         StellarModels.write_data(sm, TNUMBER)
#         StellarModels.write_terminal_info(sm)

#         update_plotter!(plotter, sm)

#         # check termination conditions
#         if (sm.props.model_number > sm.opt.termination.max_model_number)
#             StellarModels.write_terminal_info(sm; now=true)
#             println("Reached maximum model number")
#             break
#         end
#         if (exp(get_value(sm.props.lnT[1])) > sm.opt.termination.max_center_T)
#             StellarModels.write_terminal_info(sm; now=true)
#             println("Reached maximum central temperature")
#             break
#         end
#         if (get_value(sm.props.xa[1, sm.network.xa_index[:H1]]) < sm.opt.termination.min_center_H1)
#             StellarModels.write_terminal_info(sm; now=true)
#             println("Reached minimum central hydrogen concentration")
#             break
#         end

#         if (get_value(sm.props.xa[1, sm.network.xa_index[:H1]]) < sm.opt.termination.min_center_X)
#             StellarModels.write_terminal_info(sm; now=true)
#             println("Reached minimum central hydrogen abundance")
#             break
#         end

#         # get dt for coming step
#         sm.props.dt_next = get_dt_next(sm)
#     end
#     StellarModels.shut_down_IO!(sm)

#     return
# end


# """
#     get_dt_next(sm::StellarModel)

# Computes the timestep of the next evolutionary step to be taken by the StellarModel `sm` by considering all timestep
# controls (`sm.opt.timestep`).
# """
# function get_dt_next(sm::StellarModel)
#     dt_next = sm.props.dt  # this it calculated at end of step, so props.dt is the dt we used to do this step
        
#     Rsurf = exp(get_value(sm.props.lnr[sm.props.nz]))
#     Rsurf_old = exp(get_value(sm.prv_step_props.lnr[sm.prv_step_props.nz]))
#     ΔR_div_R = abs(Rsurf - Rsurf_old) / Rsurf

#     Tc = exp(get_value(sm.props.lnT[sm.props.nz]))
#     Tc_old = exp(get_value(sm.prv_step_props.lnT[sm.prv_step_props.nz]))
#     ΔTc_div_Tc = abs(Tc - Tc_old) / Tc

#     X = get_value(sm.props.xa[1, sm.network.xa_index[:H1]])
#     Xold = get_value(sm.prv_step_props.xa[1, sm.network.xa_index[:H1]])
#     ΔX = abs(X - Xold)

#     dt_nextR = dt_next * sm.opt.timestep.delta_R_limit / ΔR_div_R
#     dt_nextTc = dt_next * sm.opt.timestep.delta_Tc_limit / ΔTc_div_Tc
#     dt_nextX = dt_next * sm.opt.timestep.delta_Xc_limit / ΔX

#     min_dt = dt_next * sm.opt.timestep.dt_max_decrease
#     dt_next = min(sm.opt.timestep.dt_max_increase * dt_next, dt_nextR, dt_nextTc, dt_nextX)
#     dt_next = max(dt_next, min_dt)
#     return dt_next
# end

# """
#     compute_starting_model_properties(sm::StellarModel)

# Computes all stellar model properties based on the current state of the star.
# This sets the model up to start the evolution loop, or test the evaluation of
# the equations and the linear solver.
# """
# function compute_starting_model_properties!(sm::StellarModel)
#     StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#     StellarModels.cycle_props!(sm);
#     StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
#     StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)
#     StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
#     StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
#     StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)
# end

# """
#     do_evolution_loop(sm::StellarModel)

# Performs the main evolutionary loop of the input StellarModel `sm`. It continues taking steps until one of the
# termination criteria is reached (defined in `sm.opt.termination`).
# """
# function do_evolution_loop!(sm::StellarModel; plotter::TPLOTTER = Plotting.NullPlotter(),
#                                 max_retries_in_a_row = 10) where{TPLOTTER<:AbstractPlotter}
#     # before loop actions
#     TNUMBER = eltype(sm.props.ind_vars)  # determine type of numbers
#     StellarModels.create_output_files!(sm, TNUMBER)
#     compute_starting_model_properties!(sm)

#     # evolution loop, be sure to have sensible termination conditions or this will go on forever!
#     while true
#         StellarModels.cycle_props!(sm)  # move props of previous step to prv_step_props of current step

#         # either remesh, or copy over from prv_step_props
#         if sm.opt.remesh.do_remesh
#             StellarModels.remesher!(sm)
#         else
#             StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)
#         end
#         # time derivatives in the equations use the remeshed info, so
#         # save start_step_props before we attempt any newton solver
#         StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
#         StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
#         sm.start_step_props.dt = sm.prv_step_props.dt_next  # dt of this step becomes dt_next of previous

#         # step loop
#         sm.solver_data.newton_iters = 0
#         max_steps = sm.opt.solver.newton_max_iter
#         if (sm.start_step_props.model_number == 0)
#             max_steps = sm.opt.solver.newton_max_iter_first_step
#         end

#         StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
#         StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

#         corr = @view sm.solver_data.solver_corr[1:sm.nvars*sm.props.nz]
#         equs = @view sm.solver_data.eqs_numbers[1:sm.nvars*sm.props.nz]


#         # Exhaust timestep retries before increasing correction damping.
#         # Every trial uses the same remeshed starting model.
#         initial_dt = sm.start_step_props.dt
#         max_damping_trials = 4
#         newton_converged = false
#         for damping_attempt in 0:max_damping_trials
#             damping_factor = 0.5^damping_attempt
#             for timestep_attempt in 0:max_retries_in_a_row
#             sm.start_step_props.dt = initial_dt * sm.opt.timestep.dt_retry_decrease^timestep_attempt
#             StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
#             StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)
#             sm.solver_data.newton_iters = 0


#              # evaluate the equations for the first step
#             StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#             eval_jacobian_eqs!(sm)  # heavy lifting happens here!
#             if sm.opt.solver.report_retries
#                 @show sm.props.model_number damping_attempt damping_factor timestep_attempt sm.props.dt
#             end
#         for i = 1:max_steps
#             block_tridiagonal_solver!(sm, sm.solver_data)  # here as well

#             (abs_max_corr, i_corr) = findmax(abs, corr)
#             signed_max_corr = corr[i_corr]
#             corr_nz = div(i_corr - 1, sm.nvars) + 1
#             corr_equ = mod(i_corr - 1, sm.nvars) + 1
#             rel_corr = abs_max_corr/eps(sm.props.ind_vars[i_corr])

#             # scale correction
#             if sm.props.model_number == 0
#                 correction_limit = sm.opt.solver.initial_model_scale_max_correction
#             else
#                 correction_limit = sm.opt.solver.scale_max_correction
#             end
#             correction_multiplier = 1.0
#             for j in 1:sm.nvars
#                 sub_corr = @view sm.solver_data.solver_corr[j:sm.nvars:(sm.nvars*(sm.props.nz-1)+j)]
#                 max_sub_corr = maximum(abs, sub_corr)
#                 max_sub_corr == 0 && continue
#                 if sm.var_scaling[j] == :log || sm.var_scaling[j] == :unity
#                     correction_multiplier = min(correction_multiplier,
#                                                     correction_limit/max_sub_corr)
#                 elseif sm.var_scaling[j] == :maxval
#                     sub_ind_vars = @view sm.props.ind_vars[j:sm.nvars:(sm.nvars*(sm.props.nz-1)+j)]
#                     max_var_value = maximum(abs, sub_ind_vars)
#                     correction_multiplier = min(correction_multiplier,
#                                                     max_var_value*correction_limit/max_sub_corr)
#                 end
#             end

#             corr .*= damping_factor * correction_multiplier

#             # apply correction!
#             sm.props.ind_vars[1:sm.nvars*sm.props.nz] .+= corr[1:sm.nvars*sm.props.nz]
#             sm.solver_data.newton_iters = i

#             # evaluate the equations after correction and get residuals
#             try
#                 StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#                 eval_jacobian_eqs!(sm)  # heavy lifting happens here!

#                 (max_res, i_res) = findmax(abs, equs)
#                 res_nz = div(i_res - 1, sm.nvars) + 1
#                 res_equ = mod(i_res - 1, sm.nvars) + 1

#                 #reporting
#                 if sm.opt.solver.report_solver_progress &&
#                     i % sm.opt.solver.solver_progress_iter == 0
#                     @show sm.props.model_number, i, rel_corr, signed_max_corr, corr_nz, corr_equ, max_res, res_nz, res_equ, damping_attempt, damping_factor, timestep_attempt, sm.props.dt
#                 end
#                 #check if tolerances are satisfied
#                 if rel_corr < sm.opt.solver.relative_correction_tolerance &&
#                         max_res < sm.opt.solver.maximum_residual_tolerance
#                     if sm.props.model_number == 0
#                         println("Found first model")
#                     end
#                     newton_converged = true
#                     break  # successful, break the step loop
#                 end
#             catch e
#                 if isa(e, InterruptException)
#                     throw(e)
#                 end
#                 println("Error while evaluating equations")
#                 showerror(stdout, e)
#                 println()
#                 break
#             end
#             end  # Newton iteration loop

#             if newton_converged
#                 break
#             end
#             if sm.opt.solver.report_retries
#                 println("Newton failed with damping=$damping_factor, dt=$(sm.props.dt / SECYEAR) years")
#             end
#             end  # timestep retries

#         # Successful solve: leave the damping loop.
#         if newton_converged
#             break
#         end

#         if sm.opt.solver.report_retries
#             println("Timestep retries exhausted for damping_factor=$damping_factor")
#         end
#     end  # damping loop

#         if !newton_converged
#             StellarModels.uncycle_props!(sm)  # leave props at the last accepted model
#             println("All timestep and damping retries exhausted; terminating evolution")
#             break
#         end

#         # increment age and model number since we accept the step.
#         sm.props.time += sm.props.dt
#         sm.props.model_number += 1

#         # write state in sm.props and potential history/profiles.
#         StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#         StellarModels.write_data(sm, TNUMBER)
#         StellarModels.write_terminal_info(sm)

#         update_plotter!(plotter, sm)

#         # check termination conditions
#         if (sm.props.model_number > sm.opt.termination.max_model_number)
#             StellarModels.write_terminal_info(sm; now=true)
#             println("Reached maximum model number")
#             break
#         end
#         if (exp(get_value(sm.props.lnT[1])) > sm.opt.termination.max_center_T)
#             StellarModels.write_terminal_info(sm; now=true)
#             println("Reached maximum central temperature")
#             break
#         end

#         if (get_value(sm.props.xa[1, sm.network.xa_index[:H1]]) < sm.opt.termination.min_center_X)
#             StellarModels.write_terminal_info(sm; now=true)
#             println("Reached minimum central hydrogen abundance")
#             break
#         end

#         # get dt for coming step
#         sm.props.dt_next = get_dt_next(sm)
#     end
#     StellarModels.shut_down_IO!(sm)

#     return
# end

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


# """
#     do_evolution_loop(sm::StellarModel)

# Performs the main evolutionary loop of the input StellarModel `sm`. It continues taking steps until one of the
# termination criteria is reached (defined in `sm.opt.termination`).
# """
# function do_evolution_loop!(sm::StellarModel; plotter::TPLOTTER = Plotting.NullPlotter(), max_retries_in_a_row = 10) where {TPLOTTER<:AbstractPlotter}

#     # Before loop actions
#     TNUMBER = eltype(sm.props.ind_vars)
#     StellarModels.create_output_files!(sm, TNUMBER)
#     compute_starting_model_properties!(sm)

#     retry_count = 0
#     tdc_tolerance = 1e-5

#     # Evolution loop
#     while true

#         StellarModels.cycle_props!(sm)

#         # Temporary experiment: stop remeshing after accepted model 1240
#         allow_remesh = sm.opt.remesh.do_remesh && sm.prv_step_props.model_number < 1200

#         if allow_remesh
#             StellarModels.remesher!(sm)
#         else
#             StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)
#         end

#         # Save the previous accepted state for backward-Euler derivatives
#         StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
#         StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
#         sm.start_step_props.dt = sm.prv_step_props.dt_next

#         # Initialize the timestep attempt
#         sm.solver_data.newton_iters = 0

#         max_steps = sm.opt.solver.newton_max_iter

#         if sm.start_step_props.model_number == 0
#             max_steps = sm.opt.solver.newton_max_iter_first_step
#         end

#         retry_step = false
#         converged = false

#         StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
#         StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

#         n = sm.nvars * sm.props.nz
#         corr = @view sm.solver_data.solver_corr[1:n]
#         equs = @view sm.solver_data.eqs_numbers[1:n]

#         # Evaluate the initial trial state
#         try
#             StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#             eval_jacobian_eqs!(sm)
#         catch e
#             isa(e, InterruptException) && rethrow()
#             println("Error while evaluating initial equations")
#             showerror(stdout, e)
#             println()
#             retry_step = true
#         end

#         # Newton iterations
#         if !retry_step

#             for i in 1:max_steps

#                 try

#                     # Solve the coupled linear system
#                     block_tridiagonal_solver!(sm, sm.solver_data)

#                     # Maximum correction
#                     (abs_max_corr, i_corr) = findmax(abs, corr)
#                     signed_max_corr = corr[i_corr]

#                     corr_nz = div(i_corr-1, sm.nvars) + 1
#                     corr_equ = mod1(i_corr, sm.nvars)

#                     # Retain the existing Jems correction criterion
#                     rel_corr = abs_max_corr / eps(sm.props.ind_vars[i_corr])

#                     # Determine correction limit
#                     if sm.props.model_number == 0
#                         correction_limit = sm.opt.solver.initial_model_scale_max_correction
#                     else
#                         correction_limit = sm.opt.solver.scale_max_correction
#                     end

#                     correction_multiplier = 1.0

#                     for j in 1:sm.nvars

#                         sub_corr = @view sm.solver_data.solver_corr[j:sm.nvars:(sm.nvars*(sm.props.nz-1)+j)]
#                         max_sub_corr = maximum(abs, sub_corr)

#                         if sm.var_scaling[j] == :log || sm.var_scaling[j] == :unity

#                             correction_multiplier = min(correction_multiplier, correction_limit/max_sub_corr)

#                         elseif sm.var_scaling[j] == :maxval

#                             sub_ind_vars = @view sm.props.ind_vars[j:sm.nvars:(sm.nvars*(sm.props.nz-1)+j)]
#                             max_var_value = maximum(abs, sub_ind_vars)

#                             correction_multiplier = min(correction_multiplier, max_var_value*correction_limit/max_sub_corr)

#                         end
#                     end

#                     # Scale the correction if necessary
#                     if correction_multiplier < 1.0
#                         corr .*= correction_multiplier
#                     end

#                     # Apply Newton correction
#                     sm.props.ind_vars[1:n] .+= corr
#                     sm.solver_data.newton_iters = i

#                     # Evaluate the updated model and Jacobian
#                     StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#                     eval_jacobian_eqs!(sm)

#                     # Overall maximum residual
#                     (max_res, i_res) = findmax(abs, equs)

#                     res_nz = div(i_res-1, sm.nvars) + 1
#                     res_equ = mod1(i_res, sm.nvars)

#                     # Maximum equation-5 residual across all zones
#                     tdc_res = maximum(abs, @view equs[5:sm.nvars:n])

#                     # Solver progress
#                     if sm.opt.solver.report_solver_progress && i % sm.opt.solver.solver_progress_iter == 0
#                         @show sm.props.model_number, i, rel_corr, signed_max_corr, corr_nz, corr_equ, max_res, res_nz, res_equ, tdc_res
#                     end

#                     # Accept Newton convergence only when all criteria pass
#                     if isfinite(rel_corr) && isfinite(max_res) && isfinite(tdc_res) &&
#                         rel_corr < sm.opt.solver.relative_correction_tolerance &&
#                         max_res < sm.opt.solver.maximum_residual_tolerance &&
#                         tdc_res < tdc_tolerance

#                         converged = true

#                         if sm.props.model_number == 0
#                             println("Found first model")
#                         end

#                         break
#                     end

#                 catch e

#                     isa(e, InterruptException) && rethrow()

#                     println("Error during Newton iteration ", i)
#                     showerror(stdout, e)
#                     println()

#                     retry_step = true
#                     break

#                 end
#             end
#         end

#         # Final verification before accepting the timestep
#         if converged && !retry_step

#             try

#                 StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#                 eval_jacobian_eqs!(sm)

#                 max_res_final = maximum(abs, equs)
#                 tdc_res_final = maximum(abs, @view equs[5:sm.nvars:n])

#                 converged = isfinite(max_res_final) &&
#                     isfinite(tdc_res_final) &&
#                     max_res_final < sm.opt.solver.maximum_residual_tolerance &&
#                     tdc_res_final < tdc_tolerance

#             catch e

#                 isa(e, InterruptException) && rethrow()

#                 println("Final equation evaluation failed")
#                 showerror(stdout, e)
#                 println()

#                 converged = false

#             end
#         end

#         # Reject any timestep that has not fully converged
#         if retry_step || !converged

#             retry_count += 1

#             if retry_count > max_retries_in_a_row
#                 println("Too many retries, ending simulation")
#                 break
#             end

#             if sm.opt.solver.report_retries
#                 println("Failed to converge step $(sm.props.model_number) with timestep $(sm.props.dt/SECYEAR), retrying")
#             end

#             StellarModels.uncycle_props!(sm)
#             sm.props.dt_next *= sm.opt.timestep.dt_retry_decrease

#             continue
#         end

#         # Successful timestep
#         retry_count = 0

#         # Update age and model number
#         sm.props.time += sm.props.dt
#         sm.props.model_number += 1

#         # Write model data
#         StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#         StellarModels.write_data(sm, TNUMBER)
#         StellarModels.write_terminal_info(sm)

#         update_plotter!(plotter, sm)

#         # Termination: maximum model number
#         if sm.props.model_number > sm.opt.termination.max_model_number
#             StellarModels.write_terminal_info(sm; now=true)
#             println("Reached maximum model number")
#             break
#         end

#         # Termination: maximum central temperature
#         if exp(get_value(sm.props.lnT[1])) > sm.opt.termination.max_center_T
#             StellarModels.write_terminal_info(sm; now=true)
#             println("Reached maximum central temperature")
#             break
#         end

#         # Termination: central hydrogen depletion
#         if get_value(sm.props.xa[1, sm.network.xa_index[:H1]]) < sm.opt.termination.min_center_X
#             StellarModels.write_terminal_info(sm; now=true)
#             println("Reached minimum central hydrogen abundance")
#             break
#         end

#         # Calculate timestep for the next model
#         sm.props.dt_next = get_dt_next(sm)

#     end

#     StellarModels.shut_down_IO!(sm)

#     return
# end
# function do_evolution_loop!(sm::StellarModel; plotter::TPLOTTER = Plotting.NullPlotter(),
#                                 max_retries_in_a_row = 10) where{TPLOTTER<:AbstractPlotter}
#     # before loop actions
#     TNUMBER = eltype(sm.props.ind_vars)  # determine type of numbers
#     StellarModels.create_output_files!(sm, TNUMBER)
#     compute_starting_model_properties!(sm)
#     retry_count = 0

#     # evolution loop, be sure to have sensible termination conditions or this will go on forever!
#     while true
#         StellarModels.cycle_props!(sm)  # move props of previous step to prv_step_props of current step

#         # either remesh, or copy over from prv_step_props
#          allow_remesh = sm.opt.remesh.do_remesh && sm.prv_step_props.model_number < 1200

#         if allow_remesh
#             StellarModels.remesher!(sm)
#         else
#             StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)
#         end
#         # time derivatives in the equations use the remeshed info, so
#         # save start_step_props before we attempt any newton solver
#         StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
#         StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
#         sm.start_step_props.dt = sm.prv_step_props.dt_next  # dt of this step becomes dt_next of previous

#         # step loop
#         sm.solver_data.newton_iters = 0
#         max_steps = sm.opt.solver.newton_max_iter
#         if (sm.start_step_props.model_number == 0)
#             max_steps = sm.opt.solver.newton_max_iter_first_step
#         end

#         exit_evolution = false
#         retry_step = false
#         StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
#         StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

#         corr = @view sm.solver_data.solver_corr[1:sm.nvars*sm.props.nz]
#         equs = @view sm.solver_data.eqs_numbers[1:sm.nvars*sm.props.nz]

#         # evaluate the equations for the first step
#         StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#         eval_jacobian_eqs!(sm)  # heavy lifting happens here!
#         for i = 1:max_steps
#             block_tridiagonal_solver!(sm, sm.solver_data)  # here as well

#             (abs_max_corr, i_corr) = findmax(abs, corr)
#             signed_max_corr = corr[i_corr]
#             corr_nz = i_corr÷sm.nvars + 1
#             corr_equ = i_corr%sm.nvars
#             rel_corr = abs_max_corr/eps(sm.props.ind_vars[i_corr])

#             # scale correction
#             if sm.props.model_number == 0
#                 correction_limit = sm.opt.solver.initial_model_scale_max_correction
#             else
#                 correction_limit = sm.opt.solver.scale_max_correction
#             end
#             correction_multiplier = 1.0
#             for j in 1:sm.nvars
#                 sub_corr = @view sm.solver_data.solver_corr[j:sm.nvars:(sm.nvars*(sm.props.nz-1)+j)]
#                 max_sub_corr = maximum(abs, sub_corr)
#                 if sm.var_scaling[j] == :log || sm.var_scaling[j] == :unity
#                     correction_multiplier = min(correction_multiplier,
#                                                     correction_limit/max_sub_corr)
#                 elseif sm.var_scaling[j] == :maxval
#                     sub_ind_vars = @view sm.props.ind_vars[j:sm.nvars:(sm.nvars*(sm.props.nz-1)+j)]
#                     max_var_value = maximum(abs, sub_ind_vars)
#                     correction_multiplier = min(correction_multiplier,
#                                                     max_var_value*correction_limit/max_sub_corr)
#                 end
#             end

#             if correction_multiplier < 1.0
#                 corr .*= correction_multiplier
#             end

#             # apply correction!
#             sm.props.ind_vars[1:sm.nvars*sm.props.nz] .+= corr[1:sm.nvars*sm.props.nz]
#             sm.solver_data.newton_iters = i

#             # evaluate the equations after correction and get residuals
#             try
#                 StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#                 eval_jacobian_eqs!(sm)  # heavy lifting happens here!

#                 (max_res, i_res) = findmax(abs, equs)
#                 res_nz = i_res÷sm.nvars + 1
#                 res_equ = i_res%sm.nvars

#                 #reporting
#                 if sm.opt.solver.report_solver_progress &&
#                     i % sm.opt.solver.solver_progress_iter == 0
#                     @show sm.props.model_number, i, rel_corr, signed_max_corr, corr_nz, corr_equ, max_res, res_nz, res_equ
#                 end
#                 #check if tolerances are satisfied
#                 if rel_corr < sm.opt.solver.relative_correction_tolerance &&
#                         max_res < sm.opt.solver.maximum_residual_tolerance
#                     if sm.props.model_number == 0
#                         println("Found first model")
#                     end
#                     break  # successful, break the step loop
#                 end
#             catch e
#                 if isa(e, InterruptException)
#                     throw(e)
#                 end
#                 println("Error while evaluating equations")
#                 showerror(stdout, e)
#                 retry_step = true
#             end
#             #if not, determine if we give up or retry
#             if i == max_steps
#                 if retry_count > max_retries_in_a_row
#                     exit_evolution = true
#                     println("Too many retries, ending simulation")
#                 else
#                     retry_count = retry_count + 1
#                     retry_step = true
#                     if sm.opt.solver.report_retries
#                         println("Failed to converge step $(sm.props.model_number) with timestep $(sm.props.dt/SECYEAR), retrying")
#                     end
#                 end
#             end
#         end

#         if retry_step
#             StellarModels.uncycle_props!(sm)  # reset props to what prv_step_props contains, ie mimic state at end of previous step
#             sm.props.dt_next *= sm.opt.timestep.dt_retry_decrease # adapt dt
#             continue  # go back to top of evolution loop
#         end

#         if (exit_evolution)
#             println("Terminating evolution")
#             break
#         end

#         # step must be successful at this point
#         retry_count = 0

#         # increment age and model number since we accept the step.
#         sm.props.time += sm.props.dt
#         sm.props.model_number += 1

#         # write state in sm.props and potential history/profiles.
#         StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#         StellarModels.write_data(sm, TNUMBER)
#         StellarModels.write_terminal_info(sm)

#         update_plotter!(plotter, sm)

#         # check termination conditions
#         if (sm.props.model_number > sm.opt.termination.max_model_number)
#             StellarModels.write_terminal_info(sm; now=true)
#             println("Reached maximum model number")
#             break
#         end
#         if (exp(get_value(sm.props.lnT[1])) > sm.opt.termination.max_center_T)
#             StellarModels.write_terminal_info(sm; now=true)
#             println("Reached maximum central temperature")
#             break
#         end

#         if (get_value(sm.props.xa[1, sm.network.xa_index[:H1]]) < sm.opt.termination.min_center_X)
#             StellarModels.write_terminal_info(sm; now=true)
#             println("Reached minimum central hydrogen abundance")
#             break
#         end

#         # get dt for coming step
#         sm.props.dt_next = get_dt_next(sm)
#     end
#     StellarModels.shut_down_IO!(sm)

#     return
# end

# function do_evolution_loop!(sm::StellarModel; plotter::TPLOTTER = Plotting.NullPlotter(),
#                                 max_retries_in_a_row = 10) where{TPLOTTER<:AbstractPlotter}
#     # before loop actions
#     TNUMBER = eltype(sm.props.ind_vars)  # determine type of numbers
#     StellarModels.create_output_files!(sm, TNUMBER)
#     compute_starting_model_properties!(sm)
#     retry_count = 0 
#     correction_factors = (1.0, 0.5, 0.25, 0.2, 0.125, 0.0625)
#     correction_step = 0
#     correction_factor = correction_factors[1]
#     max_correction_steps = length(correction_factors) - 1
#     original_dt = sm.props.dt_next

#     # evolution loop, be sure to have sensible termination conditions or this will go on forever!
#     while true
#         StellarModels.cycle_props!(sm)  # move props of previous step to prv_step_props of current step #sets up sm for the current time step 

#         # either remesh, or copy over from prv_step_props
#         if sm.opt.remesh.do_remesh
#             StellarModels.remesher!(sm)
#         else
#             StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)
#         end
#         # time derivatives in the equations use the remeshed info, so
#         # save start_step_props before we attempt any newton solver
#         StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props) 
#         StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
#         sm.start_step_props.dt = sm.prv_step_props.dt_next  # dt of this step becomes dt_next of previous

#         # step loop
#         sm.solver_data.newton_iters = 0
#         max_steps = sm.opt.solver.newton_max_iter
#         if (sm.start_step_props.model_number == 0)
#             max_steps = sm.opt.solver.newton_max_iter_first_step
#         end

#         converged = false
#         retry_step = false

#         StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
#         StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

#         corr = @view sm.solver_data.solver_corr[1:sm.nvars*sm.props.nz]
#         equs = @view sm.solver_data.eqs_numbers[1:sm.nvars*sm.props.nz]

#         # evaluate the equations for the first step
#         StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#         eval_jacobian_eqs!(sm)  # heavy lifting happens here!
#         for i = 1:max_steps
#             block_tridiagonal_solver!(sm, sm.solver_data)  # here as well

#             (abs_max_corr, i_corr) = findmax(abs, corr)
#             signed_max_corr = corr[i_corr]
#             corr_nz = i_corr÷sm.nvars + 1
#             corr_equ = i_corr%sm.nvars
#             rel_corr = abs_max_corr/eps(sm.props.ind_vars[i_corr])

#             # scale correction
#             if sm.props.model_number == 0
#                 correction_limit = sm.opt.solver.initial_model_scale_max_correction
#             else
#                 correction_limit = sm.opt.solver.scale_max_correction
#             end
#             correction_limit *= correction_factor
#             correction_multiplier = 1.0
#             for j in 1:sm.nvars
#                 sub_corr = @view sm.solver_data.solver_corr[j:sm.nvars:(sm.nvars*(sm.props.nz-1)+j)]
#                 max_sub_corr = maximum(abs, sub_corr)
#                 if sm.var_scaling[j] == :log || sm.var_scaling[j] == :unity
#                     correction_multiplier = min(correction_multiplier,
#                                                     correction_limit/max_sub_corr)
#                 elseif sm.var_scaling[j] == :maxval
#                     sub_ind_vars = @view sm.props.ind_vars[j:sm.nvars:(sm.nvars*(sm.props.nz-1)+j)]
#                     max_var_value = maximum(abs, sub_ind_vars)
#                     correction_multiplier = min(correction_multiplier,
#                                                     max_var_value*correction_limit/max_sub_corr)
#                 end
#             end

#             if correction_multiplier < 1.0
#                 corr .*= correction_multiplier
#             end

#             # apply correction!
#             sm.props.ind_vars[1:sm.nvars*sm.props.nz] .+= corr[1:sm.nvars*sm.props.nz]
#             sm.solver_data.newton_iters = i

#             # evaluate the equations after correction and get residuals
#             try
#                 StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#                 eval_jacobian_eqs!(sm)  # heavy lifting happens here!

#                 (max_res, i_res) = findmax(abs, equs)
#                 res_nz = i_res÷sm.nvars + 1
#                 res_equ = i_res%sm.nvars

#                 #reporting
#                 if sm.opt.solver.report_solver_progress &&
#                     i % sm.opt.solver.solver_progress_iter == 0
#                     @show sm.props.model_number, i, rel_corr, signed_max_corr, corr_nz, corr_equ, max_res, res_nz, res_equ, retry_step, correction_step
#                 end
#                 #check if tolerances are satisfied
#                 if rel_corr < sm.opt.solver.relative_correction_tolerance &&
#                         max_res < sm.opt.solver.maximum_residual_tolerance
                    
#                     converged = true 
#                     if sm.props.model_number == 0
#                         println("Found first model")
#                     end
#                     break  # successful, break the step loop
#                 end
#             catch e
#                 if isa(e, InterruptException)
#                     throw(e)
#                 end
#                 println("Error while evaluating equations")
#                 showerror(stdout, e)
#                 retry_step = true
#                 break
#             end
#             end
#              if retry_step || !converged

#                 if retry_count < max_retries_in_a_row

#                     # Retry with a smaller timestep
#                     retry_count += 1

#                     StellarModels.uncycle_props!(sm)
#                     sm.props.dt_next *= sm.opt.timestep.dt_retry_decrease

#                     println("Timestep retry $retry_count, correction factor = $correction_factor")

#                     continue

#                 elseif correction_step < max_correction_steps

#                     # All timestep retries exhausted so halve correction limit
#                     StellarModels.uncycle_props!(sm)

#                     correction_step += 1
#                     correction_factor = correction_factors[correction_step+1]
#                     retry_count = 0

#                     sm.props.dt_next = original_dt

#                     println("Correction attempt $correction_step/$max_correction_steps, factor = $correction_factor")

#                     continue

#                 else

#                     println("All correction-limit attempts failed. Terminating evolution.")
#                     break

#                 end
#             end
#         # step must be successful at this point
#         retry_count = 0
#         correction_step = 0
#         correction_factor = 1.0 
#         # increment age and model number since we accept the step.
#         sm.props.time += sm.props.dt
#         sm.props.model_number += 1

#         # write state in sm.props and potential history/profiles.
#         StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#         StellarModels.write_data(sm, TNUMBER)
#         StellarModels.write_terminal_info(sm)

#         update_plotter!(plotter, sm)

#         # check termination conditions
#         if (sm.props.model_number > sm.opt.termination.max_model_number)
#             StellarModels.write_terminal_info(sm; now=true)
#             println("Reached maximum model number")
#             break
#         end
#         if (exp(get_value(sm.props.lnT[1])) > sm.opt.termination.max_center_T)
#             StellarModels.write_terminal_info(sm; now=true)
#             println("Reached maximum central temperature")
#             break
#         end

#         if (get_value(sm.props.xa[1, sm.network.xa_index[:H1]]) < sm.opt.termination.min_center_X)
#             StellarModels.write_terminal_info(sm; now=true)
#             println("Reached minimum central hydrogen abundance")
#             break
#         end

#         # get dt for coming step
#         sm.props.dt_next = get_dt_next(sm)
#         original_dt = sm.props.dt_next  
#     end
#     StellarModels.shut_down_IO!(sm)

#     return
# end

using LinearAlgebra

function rr_numerical_error(e)
    e isa DomainError && return true
    e isa OverflowError && return true
    e isa LinearAlgebra.SingularException && return true
    e isa LinearAlgebra.ZeroPivotException && return true
    if e isa TaskFailedException
        errors = Base.current_exceptions(e.task)
        return !isempty(errors) && all(item -> rr_numerical_error(item.exception), errors)
    elseif e isa CompositeException
        return !isempty(e.exceptions) && all(rr_numerical_error,e.exceptions)
    end
    return e isa ArgumentError && occursin("matrix contains Infs or NaNs",e.msg)
end

function rr_finite_system(sm)
    sd=sm.solver_data
    nv,nz=sm.nvars,sm.props.nz
    all(isfinite,@view(sm.props.ind_vars[1:nv*nz])) || return false
    all(isfinite,@view(sd.eqs_numbers[1:nv*nz])) || return false
    for k in 1:nz
        all(isfinite,sd.jacobian_D[k]) || return false
        k>1 && !all(isfinite,sd.jacobian_L[k]) && return false
        k<nz && !all(isfinite,sd.jacobian_U[k]) && return false
    end
    return true
end

# Detect zero rows/columns before the preconditioner divides by their maxima.
function rr_nonzero_system(sm)
    sd=sm.solver_data
    nv,nz=sm.nvars,sm.props.nz
    for k in 1:nz, a in 1:nv
        rowmax=0.0; colmax=0.0
        for b in 1:nv
            rowmax=max(rowmax,abs(sd.jacobian_D[k][a,b]))
            colmax=max(colmax,abs(sd.jacobian_D[k][b,a]))
            if k>1
                rowmax=max(rowmax,abs(sd.jacobian_L[k][a,b]))
                colmax=max(colmax,abs(sd.jacobian_U[k-1][b,a]))
            end
            if k<nz
                rowmax=max(rowmax,abs(sd.jacobian_U[k][a,b]))
                colmax=max(colmax,abs(sd.jacobian_L[k+1][b,a]))
            end
        end
        rowmax>0 && colmax>0 || return false
    end
    return true
end

# function do_evolution_loop!(sm::StellarModel;
#         plotter::TPLOTTER=Plotting.NullPlotter(),
#         max_retries_in_a_row=10,
#         remesh_pause_models=3,
#         remesh_every = 10,
#         correction_factors=(1.0,0.5,0.25,0.2,0.125,0.0625),
#         keep_successful_correction=false) where {TPLOTTER<:AbstractPlotter}

#     max_retries_in_a_row>=0 && remesh_pause_models>=0 ||
#         throw(ArgumentError("Retry and pause counts must be nonnegative"))
#     !isempty(correction_factors) && all(f->isfinite(f) && 0<f<=1,correction_factors) ||
#         throw(ArgumentError("Correction factors must lie in (0,1]"))
#     all(diff(collect(correction_factors)) .<= 0) ||
#         throw(ArgumentError("Correction factors must be nonincreasing"))
#     0<sm.opt.timestep.dt_retry_decrease<1 ||
#         throw(ArgumentError("dt_retry_decrease must lie in (0,1)"))

#     TNUMBER=eltype(sm.props.ind_vars)
#     StellarModels.create_output_files!(sm,TNUMBER)
#     trial_active=false
#     try
#         compute_starting_model_properties!(sm)
#         retry_count=0
#         correction_index=1
#         force_original_mesh=false
#         remesh_cooldown=0
#         original_dt=sm.props.dt_next

#         while true
#             StellarModels.cycle_props!(sm)
#             trial_active=true
#             attempted_dt=sm.prv_step_props.dt_next
#             isfinite(attempted_dt) && attempted_dt>0 || error("Invalid proposed timestep")
#             #tried_remesh=sm.opt.remesh.do_remesh && !force_original_mesh && remesh_cooldown==0
#             remesh_due = sm.prv_step_props.model_number % remesh_every == 0
#             tried_remesh = sm.opt.remesh.do_remesh && remesh_due && !force_original_mesh && remesh_cooldown == 0
#             correction_factor=correction_factors[correction_index]
#             converged=false
#             failure_reason=:iteration_limit

#             try
#                 if tried_remesh
#                     StellarModels.remesher!(sm)
#                 else
#                     StellarModels.copy_mesh_properties!(sm,sm.start_step_props,sm.prv_step_props)
#                 end
#                 StellarModels.copy_scalar_properties!(sm.start_step_props,sm.prv_step_props)
#                 sm.start_step_props.dt=attempted_dt
#                 StellarModels.evaluate_stellar_model_properties!(sm,sm.start_step_props)
#                 StellarModels.copy_scalar_properties!(sm.props,sm.start_step_props)
#                 StellarModels.copy_mesh_properties!(sm,sm.props,sm.start_step_props)

#                 sd=sm.solver_data
#                 sd.newton_iters=0
#                 nv=sm.nvars
#                 n=nv*sm.props.nz
#                 corr=@view sd.solver_corr[1:n]
#                 equs=@view sd.eqs_numbers[1:n]
#                 max_steps=sm.props.model_number==0 ? sm.opt.solver.newton_max_iter_first_step :
#                                                   sm.opt.solver.newton_max_iter
#                 StellarModels.evaluate_stellar_model_properties!(sm,sm.props)
#                 eval_jacobian_eqs!(sm)

#                 for i in 1:max_steps
#                     if !rr_finite_system(sm)
#                         failure_reason=:nonfinite_system
#                         break
#                     elseif !rr_nonzero_system(sm)
#                         failure_reason=:zero_jacobian_row_or_column
#                         break
#                     end
#                     block_tridiagonal_solver!(sm,sd)
#                     if !all(isfinite,corr)
#                         failure_reason=:nonfinite_correction
#                         break
#                     end
#                     abs_max_corr,i_corr=findmax(abs,corr)
#                     signed_max_corr=corr[i_corr]
#                     # Retain your original undamped ULP-based convergence check.
#                     rel_corr=abs_max_corr/eps(sm.props.ind_vars[i_corr])
#                     correction_limit=correction_factor*(sm.props.model_number==0 ?
#                         sm.opt.solver.initial_model_scale_max_correction : sm.opt.solver.scale_max_correction)
#                     isfinite(correction_limit) && correction_limit>0 ||
#                         throw(ArgumentError("Correction limit must be positive and finite"))
#                     multiplier=1.0
#                     for j in 1:nv
#                         m=maximum(abs,@view(corr[j:nv:n]))
#                         m==0 && continue
#                         if sm.var_scaling[j] == :log || sm.var_scaling[j] == :unity
#                             multiplier=min(multiplier,correction_limit/m)
#                         elseif sm.var_scaling[j] == :maxval
#                             v=maximum(abs,@view(sm.props.ind_vars[j:nv:n]))
#                             multiplier=min(multiplier,v*correction_limit/m)
#                         end
#                     end
#                     if !(isfinite(multiplier) && multiplier>0)
#                         failure_reason=:suppressed_correction
#                         break
#                     end
#                     corr .*= multiplier
#                     @views sm.props.ind_vars[1:n] .+= corr
#                     sd.newton_iters=i
#                     if !all(isfinite,@view(sm.props.ind_vars[1:n]))
#                         failure_reason=:nonfinite_state
#                         break
#                     end
#                     StellarModels.evaluate_stellar_model_properties!(sm,sm.props)
#                     eval_jacobian_eqs!(sm)
#                     if !rr_finite_system(sm)
#                         failure_reason=:nonfinite_system
#                         break
#                     end
#                     max_res,i_res=findmax(abs,equs)
#                     if sm.opt.solver.report_solver_progress && i % sm.opt.solver.solver_progress_iter==0
#                         model=sm.props.model_number
#                         corr_nz=div(i_corr-1,nv)+1; corr_equ=mod1(i_corr,nv)
#                         res_nz=div(i_res-1,nv)+1; res_equ=mod1(i_res,nv)
#                         @show model i rel_corr signed_max_corr corr_nz corr_equ max_res res_nz res_equ correction_limit tried_remesh
#                     end
#                     if rel_corr<sm.opt.solver.relative_correction_tolerance &&
#                        max_res<sm.opt.solver.maximum_residual_tolerance
#                         converged=true
#                         break
#                     end
#                 end
#             catch e
#                 # Programming errors and interrupts propagate after rollback.
#                 rr_numerical_error(e) || rethrow()
#                 failure_reason=:numerical_exception
#                 showerror(stdout,e); println()
#             end

#             if !converged
#                 # Exactly one rollback, irrespective of the chosen retry branch.
#                 StellarModels.uncycle_props!(sm)
#                 trial_active=false
#                 if tried_remesh
#                     force_original_mesh=true
#                     sm.props.dt_next=attempted_dt
#                     println("Remeshed attempt failed ($failure_reason); retry original mesh, same dt=$(attempted_dt/SECYEAR) yr, factor=$correction_factor")
#                     continue
#                 elseif retry_count<max_retries_in_a_row
#                     retry_count+=1
#                     sm.props.dt_next=attempted_dt*sm.opt.timestep.dt_retry_decrease
#                     println("Timestep retry $retry_count ($failure_reason), dt=$(sm.props.dt_next/SECYEAR) yr, factor=$correction_factor")
#                     continue
#                 elseif correction_index<length(correction_factors)
#                     correction_index+=1
#                     retry_count=0
#                     sm.props.dt_next=original_dt
#                     println("Correction retry: factor=$(correction_factors[correction_index]), dt=$(original_dt/SECYEAR) yr")
#                     continue
#                 else
#                     println("All retry attempts failed ($failure_reason). Returning last accepted model.")
#                     break
#                 end
#             end

#             if force_original_mesh
#                 remesh_cooldown=remesh_pause_models
#                 println("Original-mesh retry succeeded; skip remeshing for $remesh_pause_models accepted models")
#             elseif remesh_cooldown>0
#                 remesh_cooldown-=1
#             end
#             force_original_mesh=false
#             retry_count=0
#             if !keep_successful_correction
#                 correction_index=1
#             end
#             trial_active=false
#             sm.props.time+=sm.props.dt
#             sm.props.model_number+=1
#             # The successful Newton evaluation already refreshed properties.
#             StellarModels.write_data(sm,TNUMBER)
#             StellarModels.write_terminal_info(sm)
#             update_plotter!(plotter,sm)
#             if sm.props.model_number>sm.opt.termination.max_model_number ||
#                exp(get_value(sm.props.lnT[1]))>sm.opt.termination.max_center_T ||
#                get_value(sm.props.xa[1,sm.network.xa_index[:H1]])<sm.opt.termination.min_center_X
#                 StellarModels.write_terminal_info(sm;now=true)
#                 println("Reached an evolution termination condition")
#                 break
#             end
#             sm.props.dt_next=get_dt_next(sm)
#             original_dt=sm.props.dt_next
#         end
#     finally
#         trial_active && StellarModels.uncycle_props!(sm)
#         StellarModels.shut_down_IO!(sm)
#     end
#     return nothing
# end

function do_evolution_loop!(sm::StellarModel;
        plotter::TPLOTTER=Plotting.NullPlotter(),
        max_retries_in_a_row=10,
        remesh_pause_models=3,
        remesh_every=20,
        ms_h_depletion=0.0001,
        correction_factors=(1.0,0.5,0.25,0.2,0.125,0.0625),
        keep_successful_correction=false) where {TPLOTTER<:AbstractPlotter}

    max_retries_in_a_row >= 0 && remesh_pause_models >= 0 ||
        throw(ArgumentError("Retry and pause counts must be nonnegative"))

    remesh_every > 0 ||
        throw(ArgumentError("remesh_every must be positive"))

    isfinite(ms_h_depletion) && ms_h_depletion > 0 ||
        throw(ArgumentError("ms_h_depletion must be finite and positive"))

    !isempty(correction_factors) && all(f -> isfinite(f) && 0 < f <= 1, correction_factors) ||
        throw(ArgumentError("Correction factors must lie in (0,1]"))

    all(diff(collect(correction_factors)) .<= 0) ||
        throw(ArgumentError("Correction factors must be nonincreasing"))

    0 < sm.opt.timestep.dt_retry_decrease < 1 ||
        throw(ArgumentError("dt_retry_decrease must lie in (0,1)"))

    TNUMBER = eltype(sm.props.ind_vars)

    StellarModels.create_output_files!(sm, TNUMBER)

    trial_active = false

    try
        compute_starting_model_properties!(sm)

        retry_count = 0
        correction_index = 1

        force_original_mesh = false
        remesh_cooldown = 0

        original_dt = sm.props.dt_next

        # ---------------------------------------------------------
        # Main-sequence detection
        # ---------------------------------------------------------

        iH = sm.network.xa_index[:H1]

        Xc_initial = get_value(sm.props.xa[1, iH])

        main_sequence_started = false

        # ---------------------------------------------------------
        # Evolution loop
        # ---------------------------------------------------------

        while true

            StellarModels.cycle_props!(sm)

            trial_active = true

            attempted_dt = sm.prv_step_props.dt_next

            isfinite(attempted_dt) && attempted_dt > 0 ||
                error("Invalid proposed timestep")

            # -----------------------------------------------------
            # Decide whether to remesh
            #
            # Before MS: every remesh_every accepted models
            # During MS: every evolutionary step
            # -----------------------------------------------------

            remesh_due = main_sequence_started ||
                         sm.prv_step_props.model_number % remesh_every == 0

            tried_remesh = sm.opt.remesh.do_remesh &&
                           remesh_due &&
                           !force_original_mesh &&
                           (main_sequence_started || remesh_cooldown == 0)

            correction_factor = correction_factors[correction_index]

            converged = false
            failure_reason = :iteration_limit

            # -----------------------------------------------------
            # Attempt the stellar timestep
            # -----------------------------------------------------

            try

                if tried_remesh
                    StellarModels.remesher!(sm)
                else
                    StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)
                end

                StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)

                sm.start_step_props.dt = attempted_dt

                StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)

                StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
                StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

                sd = sm.solver_data

                sd.newton_iters = 0

                nv = sm.nvars
                n = nv * sm.props.nz

                corr = @view sd.solver_corr[1:n]
                equs = @view sd.eqs_numbers[1:n]

                max_steps = sm.props.model_number == 0 ?
                    sm.opt.solver.newton_max_iter_first_step :
                    sm.opt.solver.newton_max_iter

                StellarModels.evaluate_stellar_model_properties!(sm, sm.props)

                eval_jacobian_eqs!(sm)

                # -------------------------------------------------
                # Newton iterations
                # -------------------------------------------------

                for i in 1:max_steps

                    if !rr_finite_system(sm)

                        failure_reason = :nonfinite_system
                        break

                    elseif !rr_nonzero_system(sm)

                        failure_reason = :zero_jacobian_row_or_column
                        break

                    end

                    block_tridiagonal_solver!(sm, sd)

                    if !all(isfinite, corr)

                        failure_reason = :nonfinite_correction
                        break

                    end

                    abs_max_corr, i_corr = findmax(abs, corr)

                    signed_max_corr = corr[i_corr]

                    rel_corr = abs_max_corr / eps(sm.props.ind_vars[i_corr])

                    # ---------------------------------------------
                    # Correction limiting
                    # ---------------------------------------------

                    correction_limit = correction_factor * (
                        sm.props.model_number == 0 ?
                        sm.opt.solver.initial_model_scale_max_correction :
                        sm.opt.solver.scale_max_correction
                    )

                    isfinite(correction_limit) && correction_limit > 0 ||
                        throw(ArgumentError("Correction limit must be positive and finite"))

                    multiplier = 1.0

                    for j in 1:nv

                        m = maximum(abs, @view(corr[j:nv:n]))

                        m == 0 && continue

                        if sm.var_scaling[j] == :log || sm.var_scaling[j] == :unity

                            multiplier = min(multiplier, correction_limit/m)

                        elseif sm.var_scaling[j] == :maxval

                            v = maximum(abs, @view(sm.props.ind_vars[j:nv:n]))

                            multiplier = min(multiplier, v*correction_limit/m)

                        end

                    end

                    if !(isfinite(multiplier) && multiplier > 0)

                        failure_reason = :suppressed_correction
                        break

                    end

                    corr .*= multiplier

                    # ---------------------------------------------
                    # Apply correction
                    # ---------------------------------------------

                    @views sm.props.ind_vars[1:n] .+= corr

                    sd.newton_iters = i

                    if !all(isfinite, @view(sm.props.ind_vars[1:n]))

                        failure_reason = :nonfinite_state
                        break

                    end

                    # ---------------------------------------------
                    # Reevaluate stellar equations
                    # ---------------------------------------------

                    StellarModels.evaluate_stellar_model_properties!(sm, sm.props)

                    eval_jacobian_eqs!(sm)

                    if !rr_finite_system(sm)

                        failure_reason = :nonfinite_system
                        break

                    end

                    max_res, i_res = findmax(abs, equs)

                    # ---------------------------------------------
                    # Solver reporting
                    # ---------------------------------------------

                    if sm.opt.solver.report_solver_progress &&
                       i % sm.opt.solver.solver_progress_iter == 0

                        model = sm.props.model_number

                        corr_nz = div(i_corr-1, nv) + 1
                        corr_equ = mod1(i_corr, nv)

                        res_nz = div(i_res-1, nv) + 1
                        res_equ = mod1(i_res, nv)

                        @show model i rel_corr signed_max_corr corr_nz corr_equ max_res res_nz res_equ correction_limit tried_remesh main_sequence_started

                    end

                    # ---------------------------------------------
                    # Convergence
                    # ---------------------------------------------

                    if rel_corr < sm.opt.solver.relative_correction_tolerance &&
                       max_res < sm.opt.solver.maximum_residual_tolerance

                        converged = true
                        break

                    end

                end

            catch e

                rr_numerical_error(e) || rethrow()

                failure_reason = :numerical_exception

                showerror(stdout, e)
                println()

            end

            # -----------------------------------------------------
            # Failed timestep
            # -----------------------------------------------------

            if !converged

                StellarModels.uncycle_props!(sm)

                trial_active = false

                if tried_remesh

                    # Remeshed attempt failed.
                    # Retry original mesh at the SAME timestep.

                    force_original_mesh = true

                    sm.props.dt_next = attempted_dt

                    println("Remeshed attempt failed ($failure_reason); retry original mesh, same dt=$(attempted_dt/SECYEAR) yr, factor=$correction_factor")

                    continue

                elseif retry_count < max_retries_in_a_row

                    # Original-mesh attempt failed.
                    # Reduce timestep.

                    retry_count += 1

                    sm.props.dt_next = attempted_dt * sm.opt.timestep.dt_retry_decrease

                    println("Timestep retry $retry_count ($failure_reason), dt=$(sm.props.dt_next/SECYEAR) yr, factor=$correction_factor")

                    continue

                elseif correction_index < length(correction_factors)

                    # All timestep retries failed.
                    # Change correction factor and restore original dt.

                    correction_index += 1

                    retry_count = 0

                    sm.props.dt_next = original_dt

                    println("Correction retry: factor=$(correction_factors[correction_index]), dt=$(original_dt/SECYEAR) yr")

                    continue

                else

                    println("All retry attempts failed ($failure_reason). Returning last accepted model.")

                    break

                end

            end

            # -----------------------------------------------------
            # Successful timestep: update remeshing state
            # -----------------------------------------------------

            if force_original_mesh

                remesh_cooldown = remesh_pause_models

                if !main_sequence_started
                    println("Original-mesh retry succeeded; skip remeshing for $remesh_pause_models accepted models")
                end

            elseif remesh_cooldown > 0

                remesh_cooldown -= 1

            end

            force_original_mesh = false

            retry_count = 0

            if !keep_successful_correction
                correction_index = 1
            end

            # -----------------------------------------------------
            # Accept stellar model
            # -----------------------------------------------------

            trial_active = false

            sm.props.time += sm.props.dt
            sm.props.model_number += 1

            # -----------------------------------------------------
            # Detect beginning of main sequence
            #
            # Only check an ACCEPTED model.
            # Once enabled, the flag remains true.
            # -----------------------------------------------------

            if !main_sequence_started

                Xc = get_value(sm.props.xa[1, iH])

                if Xc_initial - Xc >= ms_h_depletion

                    main_sequence_started = true

                    remesh_cooldown = 0

                    println("Main sequence detected at model $(sm.props.model_number)")
                    println("Age = $(sm.props.time/SECYEAR) yr")
                    println("Central hydrogen: $Xc_initial → $Xc")
                    println("Enabling remeshing at every evolutionary step.")

                end

            end

            # -----------------------------------------------------
            # Write accepted model
            # -----------------------------------------------------

            StellarModels.write_data(sm, TNUMBER)

            StellarModels.write_terminal_info(sm)

            update_plotter!(plotter, sm)

            # -----------------------------------------------------
            # Termination conditions
            # -----------------------------------------------------

            if sm.props.model_number > sm.opt.termination.max_model_number ||
               exp(get_value(sm.props.lnT[1])) > sm.opt.termination.max_center_T ||
               get_value(sm.props.xa[1, iH]) < sm.opt.termination.min_center_X

                StellarModels.write_terminal_info(sm; now=true)

                println("Reached an evolution termination condition")

                break

            end

            # -----------------------------------------------------
            # Timestep for next model
            # -----------------------------------------------------

            sm.props.dt_next = get_dt_next(sm)

            original_dt = sm.props.dt_next

        end

    finally

        trial_active && StellarModels.uncycle_props!(sm)

        StellarModels.shut_down_IO!(sm)

    end

    return nothing

end