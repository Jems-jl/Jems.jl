"""
    NewtonSolver

Standard Newton solver with correction limiting. Convergence tolerances and
iteration limits are read from the model's solver options.
"""
Base.@kwdef struct NewtonSolver <: AbstractNonlinearSolver
end


"""
    solve_step!(solver, model, evolver)

Solve one implicit evolutionary step. Returns a named tuple containing
`converged`, `reason`, and `iterations`. Model-specific methods are selected by
multiple dispatch.
"""
function solve_step!(::NewtonSolver, sm::StellarModel,
                     evolver::DefaultEvolver)
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
        hook_before_evaluate_stellar_model_properties!(sm, sm.props, evolver) # hook function 
        StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
        hook_after_evaluate_stellar_model_properties!(sm, sm.props, evolver) # hook function
        eval_jacobian_eqs!(sm)  # heavy lifting happens here!]
    catch e
        e isa InterruptException && rethrow()
        println("Error while evaluating initial equations")
        showerror(stdout,e)
        println()
        return (converged=false, reason=:initial_evaluation_error, iterations=0)
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
            hook_before_evaluate_stellar_model_properties!(sm, sm.props, evolver) # hook function
            StellarModels.evaluate_stellar_model_properties!(sm, sm.props)
            hook_after_evaluate_stellar_model_properties!(sm, sm.props, evolver) # hook function
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
                return (converged=true, reason=:converged, iterations=i)
            end
        catch e
            e isa InterruptException && rethrow()
            println("Error during Newton iteration $i")
            showerror(stdout, e)
            println()
            return (converged=false, reason=:evaluation_error, iterations=i)
        end
    end
    return (converged=false, reason=:iteration_limit,
            iterations=max_steps)
end

function solve_step!(::NewtonSolver, oz::OneZone, evolver::DefaultEvolver)
    oz.solver_data.newton_iters = 0
    max_steps = oz.props.model_number == 0 ?
        oz.opt.solver.newton_max_iter_first_step :
        oz.opt.solver.newton_max_iter
    corr = oz.solver_data.solver_corr
    equs = oz.solver_data.eqs_numbers

    try
        eval_jacobian_eqs!(oz)
    catch e
        e isa InterruptException && rethrow()
        println("Error while evaluating initial equations")
        showerror(stdout, e)
        println()
        return (converged=false, reason=:initial_evaluation_error,
                iterations=0)
    end

    for i in 1:max_steps
        try
            block_tridiagonal_solver!(oz, oz.solver_data)

            abs_max_corr, i_corr = findmax(abs, corr)
            signed_max_corr = corr[i_corr]
            corr_nz = div(i_corr - 1, oz.nvars) + 1
            corr_equ = mod1(i_corr, oz.nvars)
            rel_corr = abs_max_corr / eps(oz.props.ind_vars[i_corr])

            correction_limit = oz.props.model_number == 0 ?
                oz.opt.solver.initial_model_scale_max_correction :
                oz.opt.solver.scale_max_correction
            correction_multiplier = min(1.0,
                                          correction_limit / abs_max_corr)
            correction_multiplier < 1.0 &&
                (corr .*= correction_multiplier)

            oz.props.ind_vars .+= corr
            oz.solver_data.newton_iters = i

            StellarModels.evaluate_one_zone_properties!(oz, oz.props)
            eval_jacobian_eqs!(oz)

            max_res, i_res = findmax(abs, equs)
            res_nz = div(i_res - 1, oz.nvars) + 1
            res_equ = mod1(i_res, oz.nvars)

            if oz.opt.solver.report_solver_progress &&
               i % oz.opt.solver.solver_progress_iter == 0
                @show oz.props.model_number, i, rel_corr, signed_max_corr,
                    corr_nz, corr_equ, max_res, res_nz, res_equ
            end

            if rel_corr < oz.opt.solver.relative_correction_tolerance &&
               max_res < oz.opt.solver.maximum_residual_tolerance
                oz.props.model_number == 0 && println("Found first model")
                return (converged=true, reason=:converged, iterations=i)
            end
        catch e
            e isa InterruptException && rethrow()
            println("Error during Newton iteration $i")
            showerror(stdout, e)
            println()
            return (converged=false, reason=:evaluation_error, iterations=i)
        end
    end

    return (converged=false, reason=:iteration_limit,
            iterations=max_steps)
end
