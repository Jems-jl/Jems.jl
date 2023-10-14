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
        si.lnP[i] = sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnP]]
        si.lnr[i] = sm.ind_vars[(i - 1) * sm.nvars + sm.vari[:lnr]]

        species_names = sm.var_names[(sm.nvars - sm.network.nspecies + 1):end]

        xa = view(sm.ind_vars, (i * sm.nvars - sm.network.nspecies + 1):(i * sm.nvars))

        set_EOS_resultsTP!(sm.eos, sm.psi.eos_res[i], sm.psi.lnT[i], sm.psi.lnP[i], xa, species_names)

        si.lnρ[i] = log(sm.psi.eos_res[i].ρ)
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

        dt_next = min(sm.opt.timestep.dt_max_increase * dt_next, dt_nextR, dt_nextTc, dt_nextX)
        return dt_next
    end
end

"""
    do_evolution_loop(sm::StellarModel)

Performs the main evolutionary loop of the input StellarModel `sm`. It continues taking steps until one of the
termination criteria is reached (defined in `sm.opt.termination`).
"""
function do_evolution_loop(sm::StellarModel)
    set_step_info!(sm, sm.esi)
    # evolution loop, be sure to have sensible termination conditions or this will go on forever!
    dt_factor = 1.0 # this is changed during retries to lower the timestep
    retry_count = 0
    while true
        # get dt for this step
        dt_next = get_dt_next(sm)*dt_factor

        cycle_step_info!(sm)  # move esi of previous step to psi of this step

        # remeshing
        if sm.opt.remesh.do_remesh
            sm = StellarModels.remesher!(sm)
        end

        set_step_info!(sm, sm.ssi)  # set info before we attempt any newton solver

        sm.ssi.dt = dt_next
        sm.dt = dt_next

        max_steps = sm.opt.solver.newton_max_iter
        if (sm.model_number == 0)
            max_steps = sm.opt.solver.newton_max_iter_first_step
        end

        exit_evolution = false
        retry_step = false
        # step loop
        for i = 1:max_steps
            eval_jacobian_eqs!(sm)
            thomas_algorithm!(sm)
            corr = @view sm.solver_corr[1:sm.nvars*sm.nz]

            real_max_corr = maximum(corr)

            # scale surface correction to prevent negative surface luminosity
            # if correction will produce negative L, scale it so L is halved
            corr_lum_surf = corr[sm.nvars * (sm.nz - 1) + sm.vari[:lum]]
            lum_surf = sm.ind_vars[sm.nvars * (sm.nz - 1) + sm.vari[:lum]]
            if lum_surf + corr_lum_surf < 0.0
                corr = corr * (-0.1 * lum_surf / corr_lum_surf)
            end

            # scale correction
            if sm.model_number == 0
                corr = corr * min(1, sm.opt.solver.initial_model_scale_max_correction / maximum(corr))
            else
                corr = corr * min(1, sm.opt.solver.scale_max_correction / maximum(corr))
            end
            if i % 50 == 0
                @show i, maximum(corr), real_max_corr, maximum(sm.eqs_numbers)
            end
            # first try applying correction and see if it would give negative luminosity
            for i=1:sm.nz*sm.nvars
                sm.ind_vars[i] = sm.ind_vars[i] + corr[i]
            end
            if real_max_corr < 1e-10
                if sm.model_number == 0
                    println("Found first model")
                end
                if sm.model_number % 100 == 0
                    @show sm.model_number, i, real_max_corr, maximum(sm.eqs_numbers), dt_next / SECYEAR, sm.time / SECYEAR
                end
                break
            end
            if i == max_steps
                if retry_count > 10
                    exit_evolution = true
                    println("Too many retries, ending simulation")
                else
                    retry_count = retry_count + 1
                    retry_step = true
                    println("Failed to converge step $(sm.model_number) with timestep $(dt_next/SECYEAR), retrying")
                end
            end
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
    return sm
end
