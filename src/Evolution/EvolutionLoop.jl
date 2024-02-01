"""
    set_end_step_info(sm::StellarModel)

Sets the StellarStepInfo `si`` from current state of the StellarModel `sm`.
"""
function set_step_info!(sm::StellarModel, si::StellarModels.StellarStepInfo)
    si.model_number = sm.model_number
    si.time = sm.time
    si.dt = sm.dt

    si.nz = sm.nz
    si.mstar = sm.mstar

    Threads.@threads for i = 1:(sm.nz)
        si.m[i] = sm.m[i]
        si.dm[i] = sm.dm[i]

        si.lnT[i] = sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnT]]
        si.L[i] = sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lum]]
        si.lnρ[i] = sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnρ]]
        si.lnr[i] = sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnr]]

        xa = view(sm.ind_vars, (i * sm.nvars - sm.network.nspecies + 1):(i * sm.nvars))
        si.X[i] = xa[sm.network.xa_index[:H1]]  # can later include H2 as well.
        si.Y[i] = xa[sm.network.xa_index[:He4]]  # can later include He3 as well.

        set_EOS_resultsTρ!(sm.eos, si.eos_res[i], si.lnT[i], si.lnρ[i], xa, sm.network.species_names)

        si.lnP[i] = log(si.eos_res[i].P)
        for k = 1:sm.nvars
            si.ind_vars[(i - 1) * sm.nvars + k] = sm.ind_vars[(i - 1) * sm.nvars + k]
        end
    end
end

"""
    cycle_step_info!(sm::StellarModel)

Moves the model info of the StellarModel `sm` over one state:
start step info -> end step info -> previous step info -> start step info.
"""
function cycle_step_info!(sm::StellarModel)
    temp_step_info = sm.psi
    sm.psi = sm.esi
    sm.esi = sm.ssi
    sm.ssi = temp_step_info
end

"""
    uncycle_step_info!(sm::StellarModel)

Moves the model info of the StellarModel `sm` back one state:
start step info <- end step info <- previous step info <- start step info.
"""
function uncycle_step_info!(sm::StellarModel)
    temp_step_info = sm.esi
    sm.esi = sm.psi
    sm.psi = sm.ssi
    sm.ssi = temp_step_info
end

"""
    get_dt_next(sm::StellarModel)

Computes the timestep of the next evolutionary step to be taken by the StellarModel `sm` by considering all timestep
controls (`sm.opt.timestep`).
"""
function get_dt_next(sm::StellarModel)
    dt_next = sm.esi.dt
    if (sm.esi.model_number == 0)
        return dt_next
    else
        Rsurf = exp(sm.esi.lnr[sm.esi.nz])
        Rsurf_old = exp(sm.psi.lnr[sm.esi.nz])
        ΔR_div_R = abs(Rsurf - Rsurf_old) / Rsurf

        Tc = exp(sm.esi.lnT[sm.esi.nz])
        Tc_old = exp(sm.psi.lnT[sm.esi.nz])
        ΔTc_div_Tc = abs(Tc - Tc_old) / Tc

        X = sm.esi.ind_vars[(sm.esi.nz - 1) * sm.nvars + sm.vari[:H1]]
        Xold = sm.psi.ind_vars[(sm.esi.nz - 1) * sm.nvars + sm.vari[:H1]]
        ΔX = abs(X - Xold) / (X)

        dt_nextR = dt_next * sm.opt.timestep.delta_R_limit / ΔR_div_R
        dt_nextTc = dt_next * sm.opt.timestep.delta_Tc_limit / ΔTc_div_Tc
        dt_nextX = dt_next * sm.opt.timestep.delta_Xc_limit / ΔX

        min_dt = dt_next * sm.opt.timestep.dt_max_decrease
        dt_next = min(sm.opt.timestep.dt_max_increase * dt_next, dt_nextR, dt_nextTc, dt_nextX)
        dt_next = max(dt_next, min_dt)
        return dt_next
    end
end

"""
    do_evolution_loop(sm::StellarModel)

Performs the main evolutionary loop of the input StellarModel `sm`. It continues taking steps until one of the
termination criteria is reached (defined in `sm.opt.termination`).
"""
function do_evolution_loop!(sm::StellarModel)
    set_step_info!(sm, sm.esi)
    # evolution loop, be sure to have sensible termination conditions or this will go on forever!
    dt_factor = 1.0 # this is changed during retries to lower the timestep
    retry_count = 0
    while true
        # get dt for this step
        dt_next = min(SECYEAR*1e5,get_dt_next(sm)*dt_factor)

        cycle_step_info!(sm)  # move esi of previous step to psi of this step

        # remeshing
        if sm.opt.remesh.do_remesh
            sm = StellarModels.remesher!(sm)
        end

        set_step_info!(sm, sm.ssi)  # set start step info before we attempt any newton solver

        sm.ssi.dt = dt_next
        sm.dt = dt_next
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
            corr = @view sm.solver_data.solver_corr[1:sm.nvars*sm.nz]

            ## scale surface correction to prevent negative surface luminosity
            ## if correction will produce negative L, scale it so L is halved
            #corr_lum_surf = corr[sm.nvars * (sm.nz - 1) + sm.vari[:lum]]
            #lum_surf = sm.ind_vars[sm.nvars * (sm.nz - 1) + sm.vari[:lum]]
            #if lum_surf + corr_lum_surf < 0.0
            #    corr = corr * (-0.1 * lum_surf / corr_lum_surf)
            #end

            # Find zones that have larger corrections and residual
            max_corr = -1.0 # this is the "real" correction given by the linear solver. We potentially scale it below
            max_corr_cell = -1
            max_corr_var = -1
            worst_residual = -1.0 #worst_residual_data[1]
            worst_cell = -1 #worst_residual_index÷sm.nvars + 1
            worst_equ = -1 #worst_residual_index%sm.nvars
            for k in 1:sm.nz
                for j in 1:sm.nvars
                    residual = abs(sm.solver_data.eqs_numbers[(k-1)*sm.nvars + j])
                    if residual > worst_residual
                        worst_residual = residual
                        worst_cell = k
                        worst_equ = j
                    end
                    cell_corr = abs(corr[(k-1)*sm.nvars + j])
                    if cell_corr > max_corr
                        max_corr = cell_corr
                        max_corr_cell = k
                        max_corr_var = j
                    end
                end
            end

            # scale correction
            if sm.model_number == 0
                corr = corr * min(1, sm.opt.solver.initial_model_scale_max_correction / max_corr)
            else
                corr = corr * min(1, sm.opt.solver.scale_max_correction / max_corr)
            end
            if i % 1 == 0
                @show i, sm.nz, max_corr, max_corr_cell, max_corr_var, worst_residual, worst_cell, worst_equ
            end
            # first try applying correction and see if it would give negative luminosity
            for i=1:sm.nz*sm.nvars
                sm.ind_vars[i] = sm.ind_vars[i] + corr[i]
            end
            if max_corr < 1e-5 && worst_residual < 1e-8
                if sm.model_number == 0
                    println("Found first model")
                end
                sm.solver_data.newton_iters = i
                break  # successful, break the step loop
            end
            if i == max_steps
                if retry_count > 5
                    exit_evolution = true
                    println("Too many retries, ending simulation")
                else
                    retry_count = retry_count + 1
                    #retry_step = true
                    exit_evolution = true
                    println("Failed to converge step $(sm.model_number) with timestep $(dt_next/SECYEAR), retrying")
                end
            end
            sm.solver_data.newton_iters = i
        end

        if retry_step
            dt_factor = dt_factor*0.5
            uncycle_step_info!(sm)
            continue
        else
            dt_factor = 1.0
        end

        if (exit_evolution)
            println("Terminating evolution")
            break
        end

        retry_count = 0

        # increment age and model number since we accept the step.
        sm.time = sm.time + sm.ssi.dt
        sm.model_number = sm.model_number + 1

        # write state in sm.esi and potential history/profiles.
        set_step_info!(sm, sm.esi)
        StellarModels.write_data(sm)
        StellarModels.write_terminal_info(sm)

        if sm.opt.plotting.do_plotting && sm.model_number == 1
            Plotting.init_plots!(sm)
        elseif sm.opt.plotting.do_plotting && sm.model_number % sm.opt.plotting.plotting_interval == 0
            Plotting.update_plotting!(sm)
        end

        #@show sm.model_number, sm.esi.lnP[1], sm.esi.lnP[2], sm.esi.lnP[sm.nz-1], sm.esi.lnP[sm.nz]
        #@show sm.model_number, sm.esi.lnT[1], sm.esi.lnT[2], sm.esi.lnT[sm.nz-1], sm.esi.lnT[sm.nz]
        #@show sm.dm[1], sm.dm[2], sm.dm[3]
        #@show sum(sm.dm[1:sm.nz])

        # check termination conditions
        if (sm.model_number > sm.opt.termination.max_model_number)
            println("Reached maximum model number")
            break
        end
        if (exp(sm.esi.lnT[1]) > sm.opt.termination.max_center_T)
            println("Reached maximum central temperature")
            break
        end
    end
    if sm.opt.plotting.do_plotting
        Plotting.end_of_evolution(sm)
    end

    return sm
end
