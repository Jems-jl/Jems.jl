"""
    remesher!

    Todo
"""
function remesher!(sm::StellarModel)
    do_split = Vector{Bool}(undef, sm.nz) 
    do_split .= false
    for i in 1:sm.nz-1
        a = abs(sm.ind_vars[(i-1)*sm.nvars + sm.vari[:lnP]] - sm.ind_vars[(i)*sm.nvars + sm.vari[:lnP]])
        b = sm.opt.remesh.delta_log10P_split
        if a > b
            # if the condition is satisfied, we split the largest of the two cells
            if sm.dm[i] > sm.dm[i+1]
                do_split[i] = true
            else
                do_split[i+1] = true
            end
        end
    end
    # we are ignoring the edges for now
    do_split[1] = false
    do_split[sm.nz] = false
    extra_cells = sum(do_split)
    # if allocated space is not enough, we need to reallocate everything
    if sm.nz + extra_cells > length(sm.dm)
        sm = adjusted_stellar_model_data(sm, sm.nz + extra_cells, sm.nextra);
    end
    for i=sm.nz:-1:2
        if !do_split[i]
            # move everything upwards by extra_cells
            for j in 1:sm.nvars
                sm.ind_vars[(i+extra_cells-1)*sm.nvars+j] = sm.ind_vars[(i-1)*sm.nvars+j]
            end
            sm.m[i+extra_cells] = sm.m[i]
            sm.dm[i+extra_cells] = sm.dm[i]
        else
            #split the cell and lower extra_cells by 1

            #first do cell centered values
            m_up = 0.5*(sm.m[i] + sm.m[i+1])
            if (i>2)
                m_low = 0.5*(sm.m[i-1] + sm.m[i-2])
            else
                m_low = 0.5*(sm.m[i-1])
            end
            m_new_up = 0.75*sm.m[i] + 0.25*sm.m[i-1]
            m_new_low = 0.25*sm.m[i] + 0.75*sm.m[i-1]
            m_old =  0.5*(sm.m[i] + sm.m[i-1])
            for j in 1:sm.nvars
                val_old = sm.ind_vars[(i-1)*sm.nvars+j]
                val_up = sm.ind_vars[i*sm.nvars+j]
                val_low = sm.ind_vars[(i-2)*sm.nvars+j]
                
                val_new_up = val_old + (m_new_up - m_old)/(m_up - m_old)*(val_up - val_old)
                val_new_low = val_low + (m_new_low - m_low)/(m_old - m_low)*(val_old - val_low)
                sm.ind_vars[(i+extra_cells-1)*sm.nvars+j] = val_new_up
                sm.ind_vars[(i+extra_cells-2)*sm.nvars+j] = val_new_low
            end
            sm.m[i+extra_cells] = sm.m[i]
            sm.m[i+extra_cells-1] = sm.m[i]-0.5*sm.dm[i]
            sm.dm[i+extra_cells] = 0.5*sm.dm[i]
            sm.dm[i+extra_cells-1] = 0.5*sm.dm[i]

            # need to do face centered values, right now treating everything as cell centered

            extra_cells = extra_cells - 1
        end
    end
    sm.nz = sm.nz + sum(do_split)

    return sm
end