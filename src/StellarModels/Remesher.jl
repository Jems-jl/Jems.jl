"""
    remesher!

    Todo
"""
function remesher!(sm::StellarModel)
    # we first do cell splitting
    do_split = Vector{Bool}(undef, sm.nz) 
    do_split .= false
    #maxPinstar = 0
    #maxrinstar = 0

    for i in 1:sm.nz-1
        delta_log10P = abs(sm.psi.lnP[i] - sm.psi.lnP[i+1])/log(10)
        delta_log10P_max = sm.opt.remesh.delta_log10P_max
        delta_log10r = abs(sm.psi.lnr[i] - sm.psi.lnr[i+1])/log(10)
        delta_log10r_max = sm.opt.remesh.delta_log10r_max
        delta_dm = abs(sm.psi.dm[i] - sm.psi.dm[i+1])
        #delta_dm = abs(sm.dm[i] - sm.dm[i+1])
        println(delta_dm)
        delta_dm_max = sm.opt.remesh.delta_dm_max
        if (delta_log10r > delta_log10r_max)
            # if the condition is satisfied, we split the largest of the two cells
            if sm.dm[i] > sm.dm[i+1]
                do_split[i] = true
            else
                do_split[i+1] = true
            end
        elseif (delta_log10P > delta_log10P_max) 
            # if the condition is satisfied, we split the largest of the two cells
            #println("Hey you")
            if sm.dm[i] > sm.dm[i+1]
                do_split[i] = true
            else
                do_split[i+1] = true
            end
        elseif (delta_dm > delta_dm_max) 
            # if the condition is satisfied, we split the largest of the two cells
            println("Hey hey")
            if sm.dm[i] > sm.dm[i+1]
                do_split[i] = true
            else
                do_split[i+1] = true
            end
        end
    end
    # we are ignoring the edges for now
    # do_split[sm.nz] = false
    extra_cells = sum(do_split)
    println("Extra cells = ", extra_cells)
    # if allocated space is not enough, we need to reallocate everything
    if sm.nz + extra_cells > length(sm.dm)
        println("Is this working?")
        sm = adjusted_stellar_model_data(sm, sm.nz + extra_cells, sm.nextra);
    end
    for i=sm.nz:-1:2
        if extra_cells == 0
            break
        end
        if !do_split[i]
            # move everything upwards by extra_cells
            for j in 1:sm.nvars
                sm.ind_vars[(i+extra_cells-1)*sm.nvars+j] = sm.ind_vars[(i-1)*sm.nvars+j]
            end
            # RTW: why do we need these lines here? Is mass not part of the ind_vars?
            sm.psi.m[i+extra_cells] = sm.m[i]
            sm.psi.dm[i+extra_cells] = sm.dm[i]
        else
            #split the cell and lower extra_cells by 1
            #println("Are you getting here?")
            dm_00 = sm.dm[i]
            var_00 = view(sm.ind_vars, ((i-1)*sm.nvars + 1):((i-1)*sm.nvars + sm.nvars))
            if i > 1
                dm_m1 = sm.dm[i-1]
                var_m1 = view(sm.ind_vars, ((i-2)*sm.nvars + 1):((i-2)*sm.nvars + sm.nvars))
            else
                dm_m1 = NaN
                var_m1 = []
            end
            if i < sm.nz
                dm_p1 = sm.dm[i+1]
                var_p1 = view(sm.ind_vars, ((i)*sm.nvars + 1):((i)*sm.nvars + sm.nvars))
            else
                dm_p1 = NaN
                var_p1 = []
            end
            varnew_low = view(sm.ind_vars, 
                              ((i+extra_cells-2)*sm.nvars+1):((i+extra_cells-2)*sm.nvars+sm.nvars))
            varnew_up = view(sm.ind_vars, 
                              ((i+extra_cells-1)*sm.nvars+1):((i+extra_cells-1)*sm.nvars+sm.nvars))

            #for remesh_split_function in sm.remesh_split_functions
            #    remesh_split_function(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1,
            #                            varnew_low, varnew_up)
            #end

            #println("Splitting at: ", i, " + ", extra_cells)
            #println(sm.m[i], ", ", sm.dm[i])
            #println(sm.m[i+extra_cells], ", ", sm.dm[i+extra_cells])
            sm.psi.m[i+extra_cells] = sm.m[i]
            sm.psi.m[i+extra_cells-1] = sm.m[i]-0.5*sm.dm[i]
            sm.psi.dm[i+extra_cells] = 0.5*sm.dm[i]
            sm.psi.dm[i+extra_cells-1] = 0.5*sm.dm[i]

            sm.psi.dm[i] = 0
            sm.dm[i] = 0
            #sm.dm[i] = 0 #RTW should do nothing
            #sm.m[i] = 0
            #sm.m[i+extra_cells-1] =  0
            #sm.dm[i] =  0
            #sm.psi.dm[i] = 0
            #sm.dm[i+extra_cells-1] = 0
            #println(sm.m[i+extra_cells], ", ", sm.dm[i+extra_cells])
            #println(sm.m[i], ", ", sm.dm[i])
            #println()

            extra_cells = extra_cells - 1
        end
    end
    sm.nz = sm.nz + sum(do_split)

    # we then do cell merging
    # TODO      

    return sm
end

function split_lnr_lnρ(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1,
    varnew_low, varnew_up)

    lnρ_old = var_00[sm.vari[:lnρ]] # we keep the same density on both new cells
    lnr_old = var_00[sm.vari[:lnr]]

    # radius of the lower cell is computed with the continuity equation
    # for i=1 we use the r=0 boundary condition as radius
    if i == 1
        r_face_below = 0
    else
        r_face_below = exp(var_m1[sm.vari[:lnr]])
    end
    r_low = (0.5*dm_00*3/(4π*exp(lnρ_old)) + r_face_below^3)^(1/3)

    varnew_low[sm.vari[:lnρ]] = lnρ_old
    varnew_low[sm.vari[:lnr]] = log(r_low)

    varnew_up[sm.vari[:lnρ]] = lnρ_old
    varnew_up[sm.vari[:lnr]] = lnr_old
end

function split_lum(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1,
    varnew_low, varnew_up)

    # for the center face, we use the surface boundary condition L=0
    if i==1
        L_face_below = 0
    else
        L_face_below = var_m1[sm.vari[:lum]]
    end
    L_above = var_00[sm.vari[:lum]]

    varnew_low[sm.vari[:lum]] = 0.5*(L_face_below + L_above)
    varnew_up[sm.vari[:lum]] = L_above
end

function split_lnT(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1,
    varnew_low, varnew_up)

    if i==1
        lnT_low = var_00[sm.vari[:lnT]] # Central cell remains at the center
        lnT_cell_above = var_p1[sm.vari[:lnT]]

        mnew_up = 0.75*dm_00 # mass from the core
        mcell_above = dm_00 + 0.5*dm_p1

        lnT_up = lnT_low + (lnT_cell_above - lnT_low)*mnew_up/mcell_above
    elseif i==sm.nz
        lnT_up = var_00[sm.vari[:lnT]] # Surface remians at same temperature
        lnT_cell_below = var_m1[sm.vari[:lnT]]

        mnew_low = 0.5*dm_m1 + 0.25*dm_00 # mass of new lower cell from center of lower cell
        mup = 0.5*dm_m1 + dm_00 # mass of the surface from center of lower cell

        lnT_low = lnT_cell_below + (lnT_up - lnT_cell_below)*mnew_low/mup
    else
        lnT_cell_above = var_p1[sm.vari[:lnT]]
        lnT_cell_below = var_m1[sm.vari[:lnT]]
        lnT_old = var_00[sm.vari[:lnT]]

        mnew_low = 0.5*dm_m1+0.25*dm_00 # mass at cell center of new lower cell, from center of cell below
        mold = 0.5*dm_m1+0.5*dm_00 # old mass at cell center, from center of cell below
        lnT_low = lnT_cell_below + (lnT_old - lnT_cell_below)*mnew_low/mold

        mnew_up = 0.25*dm_00 # mass from center of cell before splitting 
        mcell_above = 0.5*dm_00+0.5*dm_p1
        lnT_up = lnT_old + (lnT_cell_above - lnT_old)*mnew_up/mcell_above
    end
    varnew_low[sm.vari[:lnT]] = lnT_low
    varnew_up[sm.vari[:lnT]] = lnT_up
end

function split_xa(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1,
    varnew_low, varnew_up)
    #use same composition on both cells to preserve species
    for i in 1:sm.network.nspecies
        varnew_low[sm.nvars+1-i] = var_00[sm.nvars+1-i]
        varnew_up[sm.nvars+1-i] = var_00[sm.nvars+1-i]
    end
end