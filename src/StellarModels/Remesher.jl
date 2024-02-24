"""
    remesher!

Acts on sm.start_step_props using info from sm.prv_step_props to decide whether to merge/split cells
"""
function remesher!(sm::StellarModel)
    psp::StellarModelProperties = sm.prv_step_props  # for brevity
    ssp::StellarModelProperties = sm.start_step_props
    # we first do cell splitting
    do_split = Vector{Bool}(undef, psp.nz)
    do_split .= false
    Threads.@threads for i in 1:psp.nz-1
        a = abs(log10(get_value(psp.eos_res[i].P)) -
                log10(get_value(psp.eos_res[i+1].P)))
        b = sm.opt.remesh.delta_log10P_split
        if a > b
            # if the condition is satisfied, we split the largest of the two cells
            if psp.dm[i] > psp.dm[i+1]
                do_split[i] = true
            else
                do_split[i+1] = true
            end
        end
    end
    # we are ignoring the edges for now
    # do_split[psp.nz] = false
    extra_cells = sum(do_split)
    # if allocated space is not enough, we need to reallocate everything
    if psp.nz + extra_cells > length(psp.dm)
        adjust_props_size!(sm, psp.nz + extra_cells, sm.nextra)
    end
    for i=psp.nz:-1:1
        if !do_split[i]
            # move everything upwards by extra_cells
            for j in 1:sm.nvars
                ssp.ind_vars[(i+extra_cells-1)*sm.nvars+j] = psp.ind_vars[(i-1)*sm.nvars+j]
            end
            ssp.m[i+extra_cells] = psp.m[i]
            ssp.dm[i+extra_cells] = psp.dm[i]
        else
            # split the cell and lower extra_cells by 1
            # get old values first
            dm_00 = psp.dm[i]
            var_00 = view(psp.ind_vars, ((i-1)*sm.nvars + 1):((i-1)*sm.nvars + sm.nvars))
            if i > 1
                dm_m1 = psp.dm[i-1]
                var_m1 = view(psp.ind_vars, ((i-2)*sm.nvars + 1):((i-2)*sm.nvars + sm.nvars))
            else
                dm_m1 = NaN
                var_m1 = []
            end
            if i < psp.nz
                dm_p1 = psp.dm[i+1]
                var_p1 = view(psp.ind_vars, ((i)*sm.nvars + 1):((i)*sm.nvars + sm.nvars))
            else
                dm_p1 = NaN
                var_p1 = []
            end
            # views on the new values
            varnew_low = view(ssp.ind_vars,
                              ((i+extra_cells-2)*sm.nvars+1):((i+extra_cells-2)*sm.nvars+sm.nvars))
            varnew_up = view(ssp.ind_vars,
                              ((i+extra_cells-1)*sm.nvars+1):((i+extra_cells-1)*sm.nvars+sm.nvars))
            # populate the new values based on the old ones
            for remesh_split_function in sm.remesh_split_functions
                remesh_split_function(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
            end

            ssp.m[i+extra_cells] = psp.m[i]
            ssp.m[i+extra_cells-1] = psp.m[i]-0.5*psp.dm[i]
            ssp.dm[i+extra_cells] = 0.5*psp.dm[i]
            ssp.dm[i+extra_cells-1] = 0.5*psp.dm[i]

            extra_cells -= 1
        end
    end
    
    if (extra_cells != 0)
        throw(AssertionError("extra_cells is not zero: $(extra_cells), cell splitting has gone wrong"))
    end

    ssp.nz = psp.nz + sum(do_split)  # set the new nz

    # we then do cell merging
    # TODO

end

function split_lnr_lnρ(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)

    lnρ_old = var_00[sm.vari[:lnρ]]  # we keep the same density on both new cells
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

function split_lum(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
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

function split_lnT(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
    if i==1
        lnT_low = var_00[sm.vari[:lnT]]  # Central cell remains at the center
        lnT_cell_above = var_p1[sm.vari[:lnT]]

        mnew_up = 0.75*dm_00  # mass from the core
        mcell_above = dm_00 + 0.5*dm_p1

        lnT_up = lnT_low + (lnT_cell_above - lnT_low)*mnew_up/mcell_above
    elseif i==sm.prv_step_props.nz
        lnT_up = var_00[sm.vari[:lnT]]  # Surface remians at same temperature
        lnT_cell_below = var_m1[sm.vari[:lnT]]

        mnew_low = 0.5*dm_m1 + 0.25*dm_00  # mass of new lower cell from center of lower cell
        mup = 0.5*dm_m1 + dm_00 # mass of the surface from center of lower cell

        lnT_low = lnT_cell_below + (lnT_up - lnT_cell_below)*mnew_low/mup
    else
        lnT_cell_above = var_p1[sm.vari[:lnT]]
        lnT_cell_below = var_m1[sm.vari[:lnT]]
        lnT_old = var_00[sm.vari[:lnT]]

        mnew_low = 0.5*dm_m1+0.25*dm_00  # mass at cell center of new lower cell, from center of cell below
        mold = 0.5*dm_m1+0.5*dm_00  # old mass at cell center, from center of cell below
        lnT_low = lnT_cell_below + (lnT_old - lnT_cell_below)*mnew_low/mold

        mnew_up = 0.25*dm_00  # mass from center of cell before splitting
        mcell_above = 0.5*dm_00+0.5*dm_p1
        lnT_up = lnT_old + (lnT_cell_above - lnT_old)*mnew_up/mcell_above
    end
    varnew_low[sm.vari[:lnT]] = lnT_low
    varnew_up[sm.vari[:lnT]] = lnT_up
end

function split_xa(sm, i, dm_m1, dm_00, dm_p1, var_m1, var_00, var_p1, varnew_low, varnew_up)
    # use same composition on both cells to preserve species
    for i in 1:sm.network.nspecies
        varnew_low[sm.nvars+1-i] = var_00[sm.nvars+1-i]
        varnew_up[sm.nvars+1-i] = var_00[sm.nvars+1-i]
    end
end


"""
    adjust_props_size(sm, new_nz::Int, nextra::Int)

Returns a new StellarModelProperties object with an adjusted size. Copies over the following from the currect active
properties:

  - nz
  - dt
  - time
  - ind_vars
  - mstar
  - m
  - dm
"""
function adjust_props_size!(sm::StellarModel, new_nz::Int, nextra::Int)
    # verify that new size can contain old sm
    if sm.prv_step_props.nz > new_nz + nextra
        throw(ArgumentError("Can't fit model of size nz=$(sm.prv_step_props.nz) using new_nz=$(new_nz) and nextra=$(nextra)."))
    end
    # new properties object
    adj_props = StellarModelProperties(sm.nvars, new_nz, nextra, length(sm.network.reactions),
                                       sm.network.nspecies, sm.vari, eltype(sm.prv_step_props.ind_vars))
    # backup scalar quantities
    StellarModels.copy_scalar_properties!(adj_props, sm.prv_step_props)
    # copy the mesh qyantities (other properties are updated later)
    StellarModels.copy_mesh_properties!(sm, adj_props, sm.prv_step_props)
    sm.start_step_props = adj_props

    # also the props used later need new size:
    adj_props = StellarModelProperties(sm.nvars, new_nz, nextra, length(sm.network.reactions),
                                       sm.network.nspecies, sm.vari, eltype(sm.prv_step_props.ind_vars))
    StellarModels.copy_scalar_properties!(adj_props, sm.prv_step_props)
    StellarModels.copy_mesh_properties!(sm, adj_props, sm.prv_step_props)
    sm.props = adj_props

    # also the solver needs new arrays!
    sm.solver_data = StellarModels.SolverData(sm.nvars, new_nz, nextra, use_static_arrays,
                                              eltype(sm.prv_step_props.ind_vars))
end
