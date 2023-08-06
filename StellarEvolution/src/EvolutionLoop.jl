# set end step info (sm.esi) from current state
function set_end_step_info(sm::StellarModel)
    sm.esi.model_number = sm.model_number
    sm.esi.time = sm.time
    sm.esi.dt = sm.dt

    sm.esi.nz = sm.nz
    sm.esi.mstar = sm.mstar
    for i in 1:sm.nz
        sm.esi.m[i] = sm.m[i]
        sm.esi.dm[i] = sm.dm[i]

        sm.esi.lnT[i] = sm.ind_vars[(i-1)*sm.nvars+sm.vari[:lnT]]
        sm.esi.L[i] = sm.ind_vars[(i-1)*sm.nvars+sm.vari[:lum]]
        sm.esi.lnP[i] = sm.ind_vars[(i-1)*sm.nvars+sm.vari[:lnP]]
        sm.esi.lnr[i] = sm.ind_vars[(i-1)*sm.nvars+sm.vari[:lnr]]

        species_names = sm.varnames[sm.nvars-sm.nspecies+1:end]
        xa = sm.ind_vars[i*sm.nvars-sm.nspecies+1:i*sm.nvars]

        eos = get_EOS_resultsTP(sm.eos, sm.isotope_data, sm.psi.lnT[i], sm.psi.lnP[i], xa,species_names)
        
        sm.esi.lnρ[i] = log(eos[1])
        sm.esi.ind_vars[(i-1)*sm.nvars+1:i*sm.nvars] .= sm.ind_vars[(i-1)*sm.nvars+1:i*sm.nvars]
    end
end

function cycle_step_info(sm::StellarModel)
    temp_step_info = sm.psi
    sm.psi = sm.esi
    sm.esi = sm.ssi
    sm.ssi = temp_step_info
end

function set_start_step_info(sm::StellarModel)
    # for now, we dont do anything special before the step (ie remeshing) so we just copy things from sm.psi
    sm.ssi.model_number = sm.psi.model_number
    sm.ssi.time = sm.psi.time
    sm.ssi.dt = sm.psi.dt

    sm.ssi.nz = sm.psi.nz
    sm.ssi.mstar = sm.mstar
    for i in 1:sm.nz
        sm.ssi.m[i] = sm.psi.m[i]
        sm.ssi.dm[i] = sm.psi.dm[i]

        sm.ssi.lnT[i] = sm.psi.lnT[i]
        sm.ssi.L[i] = sm.psi.L[i]
        sm.ssi.lnP[i] = sm.psi.lnP[i]
        sm.ssi.lnr[i] = sm.psi.lnr[i]
        sm.ssi.lnρ[i] = sm.psi.lnρ[i]
        sm.ssi.ind_vars[(i-1)*sm.nvars+1:i*sm.nvars] .= sm.psi.ind_vars[(i-1)*sm.nvars+1:i*sm.nvars]
    end
end

function get_dt_next(sm::StellarModel)
    dt_next = sm.esi.dt
    if (sm.esi.model_number==0)
        return dt_next
    else
        Rsurf = exp(sm.esi.lnr[sm.esi.nz])
        Rsurf_old = exp(sm.psi.lnr[sm.esi.nz])
        ΔR_div_R = abs(Rsurf-Rsurf_old)/Rsurf
        
        Tc = exp(sm.esi.lnT[sm.esi.nz])
        Tc_old = exp(sm.psi.lnT[sm.esi.nz])
        ΔTc_div_Tc = abs(Tc-Tc_old)/Tc

        X = sm.esi.ind_vars[(sm.esi.nz-1)*sm.nvars+sm.vari[:H1]]
        Xold = sm.psi.ind_vars[(sm.esi.nz-1)*sm.nvars+sm.vari[:H1]]
        ΔX = abs(X-Xold)/(X)

        dt_nextR = dt_next*sm.opt.timestep.delta_R_limit/ΔR_div_R
        dt_nextTc = dt_next*sm.opt.timestep.delta_Tc_limit/ΔTc_div_Tc
        dt_nextX = dt_next*sm.opt.timestep.delta_Xc_limit/ΔX

        dt_next = min(2*dt_next, dt_nextR, dt_nextTc, dt_nextX)
        return dt_next
    end
end

function do_evolution_loop(sm::StellarModel)
    set_end_step_info(sm)
    #be sure to have sensible termination conditions or this will go on forever!
    while true
        dt_next = get_dt_next(sm)
    
        cycle_step_info(sm)
        set_start_step_info(sm)
    
        sm.ssi.dt = dt_next
        sm.dt = dt_next
    
        max_steps = sm.opt.solver.newton_max_iter
        if (sm.model_number==0)
            max_steps = sm.opt.solver.newton_max_iter_first_step
        end
    
        exit_evolution = false
        for i in 1:max_steps
            eval_jacobian!(sm)
            eval_eqs!(sm)
            
            sm.linear_solver.A = sm.jac
            sm.linear_solver.b = -sm.eqs
            corr =solve(sm.linear_solver)
    
            real_max_corr = maximum(corr)
            
            # scale surface correction to prevent negative surface luminosity
            # if correction will produce negative L, scale it so L is halved
            corr_lum_surf = corr[sm.nvars*(sm.nz-1)+sm.vari[:lum]]
            lum_surf = sm.ind_vars[sm.nvars*(sm.nz-1)+sm.vari[:lum]]
            if lum_surf + corr_lum_surf < 0.0
                corr = corr*(-0.1*lum_surf/corr_lum_surf)
            end
    
            #scale correction
            if sm.model_number==0
                corr = corr*min(1,3/maximum(corr))
            else
                corr = corr*min(1,0.5/maximum(corr))
            end
            if i%50==0
                @show i, maximum(corr), real_max_corr, maximum(sm.eqs)
            end
            # first try applying correction and see if it would give negative luminosity
            sm.ind_vars = sm.ind_vars+corr
            if real_max_corr<1e-10
                if sm.model_number==0
                    println("Found first model")
                end
                if sm.model_number%100==0
                    @show sm.model_number, i, real_max_corr, maximum(sm.eqs), dt_next/SECYEAR, sm.time/SECYEAR
                end
                break
            end
            if i == max_steps
                exit_evolution = true
            end
        end
    
        if (exit_evolution)
            println("Failed to converge, retry")
            break
        end
    
        sm.time = sm.time + sm.ssi.dt
        sm.model_number = sm.model_number + 1
    
        set_end_step_info(sm)
    
        write_data(sm)
    
        #check termination conditions
        if (sm.model_number > sm.opt.termination.max_model_number)
            println("Reached maximum model number")
            break
        end
        if (exp(sm.esi.lnT[1]) > sm.opt.termination.max_center_T)
            println("Reached maximum central temperature")
            break
        end
    end
end