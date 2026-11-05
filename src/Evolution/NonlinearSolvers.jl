# Nonlinear solver strategies used by the evolution drivers.
Base.@kwdef struct NewtonSolver <: AbstractNonlinearSolver
end

Base.@kwdef struct BacktrackingNewtonSolver <: AbstractNonlinearSolver
    maximum_backtracks::Int = 16
    reduction_factor::Float64 = 0.5
    armijo_constant::Float64 = 1e-4
end

StellarEvolver(; solver=NewtonSolver(), max_timestep_retries=10,
               extra_termination_conditions=()) =
    StellarEvolver(solver, max_timestep_retries, extra_termination_conditions)

OneZoneEvolver(; solver=NewtonSolver(), max_timestep_retries=10,
               extra_termination_conditions=()) =
    OneZoneEvolver(solver, max_timestep_retries, extra_termination_conditions)

function solve_step!(::NewtonSolver, sm::StellarModel, evolver::StellarEvolver,
                     state::StellarEvolutionState)
    # step loop
    sm.solver_data.newton_iters = 0
    max_steps = sm.opt.solver.newton_max_iter
    if (sm.start_step_props.model_number == 0)
        max_steps = sm.opt.solver.newton_max_iter_first_step
    end
    corr = @view sm.solver_data.solver_corr[1:sm.nvars*sm.props.nz]
    equs = @view sm.solver_data.eqs_numbers[1:sm.nvars*sm.props.nz]
    
    # evaluate the equations for the first step
    try 
        StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
        eval_jacobian_eqs!(sm)  # heavy lifting happens here!]
    catch e
        e isa InterruptException && rethrow()
        println("Error while evaluating initial equations")
        showerror(stdout,e)
        println()
        return StepResult(false, :initial_evaluation_error, 0,)
    end 

    for i = 1:max_steps
      
        try

            block_tridiagonal_solver!(sm, sm.solver_data)  # here as well

            (abs_max_corr, i_corr) = findmax(abs, corr)
            signed_max_corr = corr[i_corr]
            corr_nz = div(i_corr - 1, sm.nvars) + 1
            corr_equ = mod1(i_corr, sm.nvars)
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
            StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
            eval_jacobian_eqs!(sm)  # heavy lifting happens here!

            (max_res, i_res) = findmax(abs, equs)
            res_nz = div(i_res - 1, sm.nvars) + 1
            res_equ =mod1(i_res, sm.nvars)

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
                return StepResult(true,:converged,i)
            end
        catch e
            e isa InterruptException && rethrow()
            println("Error during Newton iteration $i")
            showerror(stdout, e)
            println()
            return StepResult(false, :evaluation_error, i)
        end
    end
    return StepResult(false, :iteration_limit, max_steps)
end

function solve_step!(solver::BacktrackingNewtonSolver, sm::StellarModel,
                     evolver::StellarEvolver, state::StellarEvolutionState)
    solver.maximum_backtracks >= 0 ||
        throw(ArgumentError("maximum_backtracks must be nonnegative"))
    0.0 < solver.reduction_factor < 1.0 ||
        throw(ArgumentError("reduction_factor must be between zero and one"))
    0.0 < solver.armijo_constant < 1.0 ||
        throw(ArgumentError("armijo_constant must be between zero and one"))

    sm.solver_data.newton_iters = 0
    max_steps = sm.start_step_props.model_number == 0 ?
        sm.opt.solver.newton_max_iter_first_step :
        sm.opt.solver.newton_max_iter
    number_of_unknowns = sm.nvars * sm.props.nz
    corr = @view sm.solver_data.solver_corr[1:number_of_unknowns]
    equs = @view sm.solver_data.eqs_numbers[1:number_of_unknowns]
    ind_vars = @view sm.props.ind_vars[1:number_of_unknowns]

    try
        StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
        eval_jacobian_eqs!(sm)
    catch e
        e isa InterruptException && rethrow()
        println("Error while evaluating initial equations")
        showerror(stdout, e)
        println()
        return StepResult(false, :initial_evaluation_error, 0)
    end

    for i in 1:max_steps
        variables_before = copy(ind_vars)
        merit_before = 0.5 * sum(abs2, equs)

        try
            block_tridiagonal_solver!(sm, sm.solver_data)
        catch e
            e isa InterruptException && rethrow()
            println("Error while solving the Newton system at iteration $i")
            showerror(stdout, e)
            println()
            return StepResult(false, :linear_solve_error, i)
        end

        correction_limit = sm.props.model_number == 0 ?
            sm.opt.solver.initial_model_scale_max_correction :
            sm.opt.solver.scale_max_correction
        correction_multiplier = 1.0
        for j in 1:sm.nvars
            sub_corr = @view corr[j:sm.nvars:number_of_unknowns]
            max_sub_corr = maximum(abs, sub_corr)
            if sm.var_scaling[j] == :log || sm.var_scaling[j] == :unity
                correction_multiplier = min(
                    correction_multiplier, correction_limit / max_sub_corr)
            elseif sm.var_scaling[j] == :maxval
                sub_ind_vars = @view ind_vars[j:sm.nvars:number_of_unknowns]
                max_var_value = maximum(abs, sub_ind_vars)
                correction_multiplier = min(
                    correction_multiplier,
                    max_var_value * correction_limit / max_sub_corr)
            end
        end
        correction_multiplier < 1.0 && (corr .*= correction_multiplier)

        # Use the limited Newton direction for the convergence test. Testing the
        # damped correction could falsely report convergence merely because the
        # line search selected a very small step length.
        abs_max_corr, i_corr = findmax(abs, corr)
        signed_max_corr = corr[i_corr]
        corr_nz = div(i_corr - 1, sm.nvars) + 1
        corr_equ = mod1(i_corr, sm.nvars)
        rel_corr = abs_max_corr / abs(eps(sm.props.ind_vars[i_corr]))

        accepted = false
        accepted_alpha = 0.0
        max_res = Inf
        i_res = 1
        alpha = 1.0

        for _ in 0:solver.maximum_backtracks
            ind_vars .= variables_before .+ alpha .* corr
            try
                StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
                eval_jacobian_eqs!(sm)
                trial_merit = 0.5 * sum(abs2, equs)
                trial_max_res, trial_i_res = findmax(abs, equs)
                sufficient_decrease = (1.0 - solver.armijo_constant * alpha) * merit_before
                residual_converged = trial_max_res <
                    sm.opt.solver.maximum_residual_tolerance
                if isfinite(trial_merit) &&
                   (residual_converged || trial_merit <= sufficient_decrease)
                    accepted = true
                    accepted_alpha = alpha
                    max_res = trial_max_res
                    i_res = trial_i_res
                    break
                end
            catch e
                e isa InterruptException && rethrow()
                # An invalid trial state is treated like a rejected line-search
                # point. The next trial starts again from variables_before.
            end
            alpha *= solver.reduction_factor
        end

        if !accepted
            ind_vars .= variables_before
            return StepResult(false, :line_search_failed, i)
        end

        sm.solver_data.newton_iters = i
        res_nz = div(i_res - 1, sm.nvars) + 1
        res_equ = mod1(i_res, sm.nvars)

        if sm.opt.solver.report_solver_progress &&
           i % sm.opt.solver.solver_progress_iter == 0
            @show sm.props.model_number, i, rel_corr, signed_max_corr,
                corr_nz, corr_equ, max_res, res_nz, res_equ, accepted_alpha
        end

        if rel_corr < sm.opt.solver.relative_correction_tolerance &&
           max_res < sm.opt.solver.maximum_residual_tolerance
            sm.props.model_number == 0 && println("Found first model")
            return StepResult(true, :converged, i)
        end
    end

    return StepResult(false, :iteration_limit, max_steps)
end
