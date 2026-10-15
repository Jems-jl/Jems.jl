

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
    Tc_old = exp(get_value(sm.prv_step_props.lnT[sm.props.nz]))
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


"""
Stellar evolution loop for basic mlt 
"""
function do_evolution_loop!(sm::StellarModel; plotter::TPLOTTER = Plotting.NullPlotter(),
                                max_retries_in_a_row = 10) where{TPLOTTER<:AbstractPlotter}
    # before loop actions
    TNUMBER = eltype(sm.props.ind_vars)  # determine type of numbers
    StellarModels.create_output_files!(sm, TNUMBER)
    compute_starting_model_properties!(sm)
    retry_count = 0

    # evolution loop, be sure to have sensible termination conditions or this will go on forever!
    while true
        StellarModels.cycle_props!(sm)  # move props of previous step to prv_step_props of current step

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
            block_tridiagonal_solver!(sm, sm.solver_data)  # here as well

            (abs_max_corr, i_corr) = findmax(abs, corr)
            signed_max_corr = corr[i_corr]
            corr_nz = i_corr÷sm.nvars + 1
            corr_equ = i_corr%sm.nvars
            rel_corr = abs_max_corr/eps(sm.props.ind_vars[i_corr])

            # scale correction
            if sm.props.model_number == 0
                correction_limit = sm.opt.solver.initial_model_scale_max_correction
            else
                correction_limit = sm.opt.solver.scale_max_correction
            end
            correction_multiplier = 1.0
            for j in 1:sm.nvars
                sub_corr = @view sm.solver_data.solver_corr[j:sm.nvars:(sm.nvars*(sm.props.nz-1)+j)]
                max_sub_corr = maximum(abs, sub_corr)
                if sm.var_scaling[j] == :log || sm.var_scaling[j] == :unity
                    correction_multiplier = min(correction_multiplier,
                                                    correction_limit/max_sub_corr)
                elseif sm.var_scaling[j] == :maxval
                    sub_ind_vars = @view sm.props.ind_vars[j:sm.nvars:(sm.nvars*(sm.props.nz-1)+j)]
                    max_var_value = maximum(abs, sub_ind_vars)
                    correction_multiplier = min(correction_multiplier,
                                                    max_var_value*correction_limit/max_sub_corr)
                end
            end

            if correction_multiplier < 1.0
                corr .*= correction_multiplier
            end

            # apply correction!
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
                if retry_count > max_retries_in_a_row
                    exit_evolution = true
                    println("Too many retries, ending simulation")
                else
                    retry_count = retry_count + 1
                    retry_step = true
                    if sm.opt.solver.report_retries
                        println("Failed to converge step $(sm.props.model_number) with timestep $(sm.props.dt/SECYEAR), retrying")
                    end
                end
            end
        end

        if retry_step
            StellarModels.uncycle_props!(sm)  # reset props to what prv_step_props contains, ie mimic state at end of previous step
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
        StellarModels.write_data(sm, TNUMBER)
        StellarModels.write_terminal_info(sm)

        update_plotter!(plotter, sm)

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

        if (get_value(sm.props.xa[1, sm.network.xa_index[:H1]]) < sm.opt.termination.min_center_X)
            StellarModels.write_terminal_info(sm; now=true)
            println("Reached minimum central hydrogen abundance")
            break
        end

        # get dt for coming step
        sm.props.dt_next = get_dt_next(sm)
    end
    StellarModels.shut_down_IO!(sm)

    return
end

"""
Stellar evolution loop for step overshooting 
"""
# @inline function build_step_overshoot!(
#     props;
#     alpha_ov,
#     f0 = 0.005
# )

#     fill!(props.D_overshoot, 0.0)

#     nzf = props.nz - 1
#     nzf < 3 && return

#     # ============================================================
#     # 1. Find formal Schwarzschild boundary
#     # ============================================================

#     discriminant(i) = get_value(props.turb_res[i].∇ᵣ) - get_value(props.∇ₐ_face[i])
#     # Convective core must begin at centre.
#     d_prev = discriminant(1)
#     if !isfinite(d_prev) || d_prev <= 0
#         return
#     end

#     i_rad = nothing
#     for i in 2:nzf
#         d = discriminant(i)
#         isfinite(d) || return
#         if d <= 0
#             i_rad = i
#             break
#         end
#         d_prev = d
#     end
#     isnothing(i_rad) && return

#     i_conv = i_rad - 1

#     d_conv = discriminant(i_conv)
#     d_rad  = discriminant(i_rad)

#     # Fraction between convective and radiative face where
#     # nabla_rad - nabla_ad = 0.
#     w_b = d_conv / (d_conv - d_rad)
#     if !isfinite(w_b) || !(0.0 <= w_b <= 1.0)
#         return
#     end
#     # ============================================================
#     # 2. Interpolate physical quantities to Schwarzschild boundary
#     # ============================================================

#     logr_conv = get_value(props.lnr[i_conv])
#     logr_rad  = get_value(props.lnr[i_rad])

#     logP_conv = get_value(props.lnP_face[i_conv])
#     logP_rad  = get_value(props.lnP_face[i_rad])

#     logρ_conv = get_value(props.lnρ_face[i_conv])
#     logρ_rad  = get_value(props.lnρ_face[i_rad])

#     m_conv = props.m[i_conv]
#     m_rad  = props.m[i_rad]

#     logr_b = logr_conv + w_b * (logr_rad - logr_conv)
#     logP_b = logP_conv + w_b * (logP_rad - logP_conv)
#     logρ_b = logρ_conv + w_b * (logρ_rad - logρ_conv)
#     m_b    = m_conv    + w_b * (m_rad - m_conv)

#     r_b = exp(logr_b)
#     P_b = exp(logP_b)
#     ρ_b = exp(logρ_b)

#     if !all(isfinite, (r_b, P_b, ρ_b, m_b)) ||
#        r_b <= 0 || P_b <= 0 || ρ_b <= 0 || m_b <= 0
#         return
#     end

#     # ============================================================
#     # 3. Hp at formal convective boundary
#     # ============================================================

#     g_b = CGRAV * m_b / r_b^2
#     Hp_b = P_b / (ρ_b * g_b)

#     if !isfinite(Hp_b) || Hp_b <= 0
#         return
#     end
#     # Geometric limiter
#     H_eff = min(Hp_b, r_b)

#     # ============================================================
#     # 4. MESA-style geometry
#     #
#     # r0: point f0 Hp inside the formal boundary
#     #
#     # f_total = alpha_ov + f0
#     #
#     # r_out = r0 + f_total Hp
#     #
#     # therefore:
#     #
#     # r_out - r_b = alpha_ov Hp
#     # ============================================================

#     f_total = alpha_ov + f0

#     r0 = r_b - f0 * H_eff
#     r_out = r0 + f_total * H_eff

#     r0 = max(r0, 0.0)

#     # ============================================================
#     # 5 & 6. Find bracket and interpolate D_MLT at r0
#     # ============================================================
#     r_face_1 = exp(get_value(props.lnr[1]))
#     D0 = 0.0

#     if r0 <= r_face_1
#         # r0 is inside the innermost cell. No interpolation needed, 
#         # just use the innermost diffusion coefficient.
#         D0 = get_value(props.D_turb[1])
#     else
#         i0 = nothing
#         for i in 1:(i_conv)
#             r_left  = exp(get_value(props.lnr[i]))
#             r_right = exp(get_value(props.lnr[i+1]))

#             if r_left <= r0 <= r_right
#                 i0 = i
#                 break
#             end
#         end

#         isnothing(i0) && return

#         r_left  = exp(get_value(props.lnr[i0]))
#         r_right = exp(get_value(props.lnr[i0+1]))
#         D_left  = get_value(props.D_turb[i0])
#         D_right = get_value(props.D_turb[i0+1])

#        if D_left > 0 && D_right > 0
#             w0 = (r0 - r_left) / (r_right - r_left)
#             D0 = exp(log(D_left) + w0 * (log(D_right) - log(D_left)))
#         elseif D_left > 0
#             D0 = D_left
#         else
#             return
#         end
#     end

#     if !isfinite(D0) || D0 <= 0
#         return
#     end

#     # ============================================================
#     # 7. Apply the constant step FROM r0, not from rb
#     # ============================================================

#     for i in 1:nzf

#         r_face = exp(get_value(props.lnr[i]))

#         if r0 <= r_face <= r_out

#             props.D_overshoot[i] = D0

#         elseif r_face > r_out

#             break

#         end
#     end

#     return
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
#         fill!(sm.start_step_props.D_overshoot, 0.0)
#         sm.start_step_props.dt = sm.prv_step_props.dt_next  # dt of this step becomes dt_next of previous
#         StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
#         build_step_overshoot!(sm.start_step_props; alpha_ov = 0.0, f0 = 0.0)
#         StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
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

#         iH = sm.network.xa_index[:H1]

#         Xc_new = get_value(sm.props.xa[1, iH])
#         Xc_old = get_value(sm.start_step_props.xa[1, iH])
#         delta_Xc = abs(Xc_new - Xc_old)

#         if delta_Xc > sm.opt.timestep.delta_Xc_limit
#             retry_step = true
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


"""
Stellar evolution loop for TDC
"""
# using LinearAlgebra

# function rr_numerical_error(e)
#     e isa DomainError && return true
#     e isa OverflowError && return true
#     e isa LinearAlgebra.SingularException && return true
#     e isa LinearAlgebra.ZeroPivotException && return true
#     if e isa TaskFailedException
#         errors = Base.current_exceptions(e.task)
#         return !isempty(errors) && all(item -> rr_numerical_error(item.exception), errors)
#     elseif e isa CompositeException
#         return !isempty(e.exceptions) && all(rr_numerical_error,e.exceptions)
#     end
#     return e isa ArgumentError && occursin("matrix contains Infs or NaNs",e.msg)
# end

# function rr_finite_system(sm)
#     sd=sm.solver_data
#     nv,nz=sm.nvars,sm.props.nz
#     all(isfinite,@view(sm.props.ind_vars[1:nv*nz])) || return false
#     all(isfinite,@view(sd.eqs_numbers[1:nv*nz])) || return false
#     for k in 1:nz
#         all(isfinite,sd.jacobian_D[k]) || return false
#         k>1 && !all(isfinite,sd.jacobian_L[k]) && return false
#         k<nz && !all(isfinite,sd.jacobian_U[k]) && return false
#     end
#     return true
# end

# # Detect zero rows/columns before the preconditioner divides by their maxima.
# function rr_nonzero_system(sm)
#     sd=sm.solver_data
#     nv,nz=sm.nvars,sm.props.nz
#     for k in 1:nz, a in 1:nv
#         rowmax=0.0; colmax=0.0
#         for b in 1:nv
#             rowmax=max(rowmax,abs(sd.jacobian_D[k][a,b]))
#             colmax=max(colmax,abs(sd.jacobian_D[k][b,a]))
#             if k>1
#                 rowmax=max(rowmax,abs(sd.jacobian_L[k][a,b]))
#                 colmax=max(colmax,abs(sd.jacobian_U[k-1][b,a]))
#             end
#             if k<nz
#                 rowmax=max(rowmax,abs(sd.jacobian_U[k][a,b]))
#                 colmax=max(colmax,abs(sd.jacobian_L[k+1][b,a]))
#             end
#         end
#         rowmax>0 && colmax>0 || return false
#     end
#     return true
# end
# function do_evolution_loop!(sm::StellarModel;
#         plotter::TPLOTTER=Plotting.NullPlotter(), max_retries_in_a_row=10,
#         remesh_pause_models=3, remesh_every=5,
#         correction_factors=(1.0, 0.5, 0.25, 0.2, 0.125, 0.0625),
#         keep_successful_correction=false) where {TPLOTTER<:AbstractPlotter}

#     (max_retries_in_a_row >= 0 && remesh_pause_models >= 0) || throw(ArgumentError("Retry/pause counts must be >= 0"))
#     remesh_every > 0 || throw(ArgumentError("remesh_every must be > 0"))
#     (!isempty(correction_factors) && all(f -> isfinite(f) && 0 < f <= 1, correction_factors)) || throw(ArgumentError("Correction factors must lie in (0,1]"))
#     all(diff(collect(correction_factors)) .<= 0) || throw(ArgumentError("Correction factors must be nonincreasing"))
#     (0 < sm.opt.timestep.dt_retry_decrease < 1) || throw(ArgumentError("dt_retry_decrease must lie in (0,1)"))

#     TNUMBER = eltype(sm.props.ind_vars)
#     StellarModels.create_output_files!(sm, TNUMBER)
#     trial_active = false

#     try
#         compute_starting_model_properties!(sm)
#         retry_count, correction_index, remesh_cooldown = 0, 1, 0
#         force_original_mesh = false
#         original_dt = sm.props.dt_next
#         iH = sm.network.xa_index[:H1]

#         while true
#             StellarModels.cycle_props!(sm)
#             trial_active = true
#             attempted_dt = sm.prv_step_props.dt_next
#             (isfinite(attempted_dt) && attempted_dt > 0) || error("Invalid proposed timestep")

#             remesh_due = sm.prv_step_props.model_number % remesh_every == 0
#             tried_remesh = sm.opt.remesh.do_remesh && remesh_due && !force_original_mesh && remesh_cooldown == 0
#             correction_factor = correction_factors[correction_index]
#             converged, failure_reason = false, :iteration_limit

#             try
#                 tried_remesh ? StellarModels.remesher!(sm) : StellarModels.copy_mesh_properties!(sm, sm.start_step_props, sm.prv_step_props)
#                 StellarModels.copy_scalar_properties!(sm.start_step_props, sm.prv_step_props)
#                 sm.start_step_props.dt = attempted_dt

#                 StellarModels.evaluate_stellar_model_properties!(sm, sm.start_step_props)
#                 StellarModels.copy_scalar_properties!(sm.props, sm.start_step_props)
#                 StellarModels.copy_mesh_properties!(sm, sm.props, sm.start_step_props)

#                 sd, nv = sm.solver_data, sm.nvars
#                 n = nv * sm.props.nz
#                 sd.newton_iters = 0
#                 corr, equs = @view(sd.solver_corr[1:n]), @view(sd.eqs_numbers[1:n])
#                 max_steps = sm.props.model_number == 0 ? sm.opt.solver.newton_max_iter_first_step : sm.opt.solver.newton_max_iter

#                 StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#                 eval_jacobian_eqs!(sm)

#                 for i in 1:max_steps
#                     rr_finite_system(sm) || (failure_reason = :nonfinite_system; break)
#                     rr_nonzero_system(sm) || (failure_reason = :zero_jacobian_row_or_column; break)

#                     block_tridiagonal_solver!(sm, sd)
#                     all(isfinite, corr) || (failure_reason = :nonfinite_correction; break)

#                     abs_max_corr, i_corr = findmax(abs, corr)
#                     rel_corr = abs_max_corr / eps(sm.props.ind_vars[i_corr])

#                     limit_base = sm.props.model_number == 0 ? sm.opt.solver.initial_model_scale_max_correction : sm.opt.solver.scale_max_correction
#                     correction_limit = correction_factor * limit_base
#                     (isfinite(correction_limit) && correction_limit > 0) || throw(ArgumentError("Correction limit must be positive and finite"))

#                     multiplier = 1.0
#                     for j in 1:nv
#                         m = maximum(abs, @view(corr[j:nv:n]))
#                         m == 0 && continue
#                         s_type = sm.var_scaling[j]
                        
#                         if s_type ∈ (:log, :unity)
#                             multiplier = min(multiplier, correction_limit / m)
#                         elseif s_type == :maxval
#                             v = maximum(abs, @view(sm.props.ind_vars[j:nv:n]))
#                             multiplier = min(multiplier, v * correction_limit / m)
#                         end
#                     end

#                     # Abundance boundary limiter
#                     nspecies, nhydro = sm.network.nspecies, nv - sm.network.nspecies
#                     boundary_fraction, X_abs_tol = 0.9, 1e-14
#                     abundance_multiplier, abundance_failure = 1.0, false
#                     limiting_zone, limiting_species, limiting_reason = 0, 0, :none

#                     for k in 1:sm.props.nz
#                         for j in 1:nspecies
#                             idx = (k - 1) * nv + nhydro + j
#                             X, dX = sm.props.ind_vars[idx], corr[idx]

#                             if !isfinite(X) || !isfinite(dX) || X < 0.0 || X > 1.0
#                                 failure_reason, abundance_failure = :invalid_abundance, true
#                                 break
#                             end

#                             Xtrial = X + dX
#                             0.0 <= Xtrial <= 1.0 && continue

#                             if Xtrial < 0.0
#                                 if X <= X_abs_tol && Xtrial >= -X_abs_tol
#                                     corr[idx] = -X
#                                     continue
#                                 end
#                                 if dX < 0.0
#                                     allowed_multiplier = boundary_fraction * X / (-dX)
#                                     if allowed_multiplier < abundance_multiplier
#                                         abundance_multiplier, limiting_zone, limiting_species, limiting_reason = allowed_multiplier, k, j, :lower_boundary
#                                     end
#                                 else
#                                     failure_reason, abundance_failure = :inconsistent_abundance_correction, true
#                                     break
#                                 end
#                             elseif Xtrial > 1.0
#                                 if dX > 0.0
#                                     allowed_multiplier = boundary_fraction * (1.0 - X) / dX
#                                     if allowed_multiplier < abundance_multiplier
#                                         abundance_multiplier, limiting_zone, limiting_species, limiting_reason = allowed_multiplier, k, j, :upper_boundary
#                                     end
#                                 else
#                                     failure_reason, abundance_failure = :inconsistent_abundance_correction, true
#                                     break
#                                 end
#                             end
#                         end
#                         abundance_failure && break
#                     end
#                     abundance_failure && break

#                     if abundance_multiplier < multiplier && sm.opt.solver.report_solver_progress
#                         @show i, limiting_zone, sm.network.species_names[limiting_species], limiting_reason, abundance_multiplier
#                     end

#                     multiplier = min(multiplier, abundance_multiplier)
#                     (isfinite(multiplier) && multiplier > 0.0) || (failure_reason = :suppressed_correction; break)

#                     corr .*= multiplier
#                     @views sm.props.ind_vars[1:n] .+= corr
#                     sd.newton_iters = i

#                     all(isfinite, @view(sm.props.ind_vars[1:n])) || (failure_reason = :nonfinite_state; break)

#                     StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
#                     eval_jacobian_eqs!(sm)
#                     rr_finite_system(sm) || (failure_reason = :nonfinite_system; break)

#                     max_res, i_res = findmax(abs, equs)

#                     if sm.opt.solver.report_solver_progress && i % sm.opt.solver.solver_progress_iter == 0
#                         @show sm.props.model_number, i, rel_corr, corr[i_corr], div(i_corr-1, nv)+1, mod1(i_corr, nv), max_res, div(i_res-1, nv)+1, mod1(i_res, nv), correction_limit, tried_remesh
#                     end

#                     if rel_corr < sm.opt.solver.relative_correction_tolerance && max_res < sm.opt.solver.maximum_residual_tolerance
#                         converged = true
#                         break
#                     end
#                 end
#             catch e
#                 rr_numerical_error(e) || rethrow()
#                 failure_reason = :numerical_exception
#                 showerror(stdout, e)
#                 println()
#             end

#             if !converged
#                 StellarModels.uncycle_props!(sm)
#                 trial_active = false

#                 if tried_remesh
#                     force_original_mesh = true
#                     sm.props.dt_next = attempted_dt
#                     println("Remeshed attempt failed ($failure_reason); retry original mesh, same dt=$(attempted_dt/SECYEAR) yr, factor=$correction_factor")
#                     continue
#                 elseif retry_count < max_retries_in_a_row
#                     retry_count += 1
#                     sm.props.dt_next = attempted_dt * sm.opt.timestep.dt_retry_decrease
#                     println("Timestep retry $retry_count ($failure_reason), dt=$(sm.props.dt_next/SECYEAR) yr, factor=$correction_factor")
#                     continue
#                 elseif correction_index < length(correction_factors)
#                     correction_index += 1
#                     retry_count = 0
#                     sm.props.dt_next = original_dt
#                     println("Correction retry: factor=$(correction_factors[correction_index]), dt=$(original_dt/SECYEAR) yr")
#                     continue
#                 else
#                     println("All retry attempts failed ($failure_reason). Returning last accepted model.")
#                     break
#                 end
#             end

#             if force_original_mesh
#                 remesh_cooldown = remesh_pause_models
#                 println("Original-mesh retry succeeded; skip remeshing for $remesh_pause_models accepted models")
#             elseif remesh_cooldown > 0
#                 remesh_cooldown -= 1
#             end

#             force_original_mesh, retry_count = false, 0
#             !keep_successful_correction && (correction_index = 1)
#             trial_active = false

#             sm.props.time += sm.props.dt
#             sm.props.model_number += 1

#             StellarModels.write_data(sm, TNUMBER)
#             StellarModels.write_terminal_info(sm)
#             update_plotter!(plotter, sm)

#             term = sm.opt.termination
#             if sm.props.model_number > term.max_model_number || exp(get_value(sm.props.lnT[1])) > term.max_center_T || get_value(sm.props.xa[1, iH]) < term.min_center_X
#                 StellarModels.write_terminal_info(sm; now=true)
#                 println("Reached an evolution termination condition")
#                 break
#             end

#             sm.props.dt_next = get_dt_next(sm)
#             original_dt = sm.props.dt_next
#         end
#     finally
#         trial_active && StellarModels.uncycle_props!(sm)
#         StellarModels.shut_down_IO!(sm)
#     end
#     return nothing
# end
